#################################################################
##                         SMART_FPinder                      ###
##                Extremely experimental version              ###
##   Developed by Hyunwoo Kim, Chen Zhang, and Raphael Reher  ###
##         William H. Gerwich and Garrison W. Cottrell        ###
##                        October, 2019                       ###
#################################################################


import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import datasets, layers, models

import numpy as np
import pandas as pd 
import datetime
import os
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs

import matplotlib as mpl
import matplotlib.pyplot as plt
from math import sqrt
import time
import json

import argparse
scale = 128 #scale of the input data. This scale will not be changed.

def load_models(models_folder="models"):
    #importing trained model
    #If no gpus are available, these models are working with CPU automatically
    #each channel has two models
    with tf.device('/CPU:0'):
        model_1ch = keras.models.load_model(os.path.join(models_folder, '(012521)SMART3_v3_1ch_multitask.hdf5')) #To predict chemical fingerprints and molecular weights
        model_2ch = keras.models.load_model(os.path.join(models_folder, '(012521)SMART3_v3_2ch_multitask.hdf5')) #To predict chemical fingerprints and molecular weights
        
    return model_1ch, model_2ch 

def load_db(db_folder="."):
    #loading DB
    
    DB =np.array(pd.read_json('DB_010621_SM3.json'))
    with open('superclass.json','r') as r:
        index_super = json.load(r)
    return DB, index_super

#This is a binary cosine between two sets
def cosine(x,y):
    '''x, y are same shape array'''
    a = set(x)
    b = set(y)
    return (len(a&b))/(sqrt(len(a))*sqrt(len(b)))

# This function is for converting CSV, TSV, or xlsx file to numpy array (128 x 128 x 1) or (128 x 128 x 2)
def CSV_converter(CSV_converter, channel): 
    _, ext = os.path.splitext(CSV_converter) # extracting filename extension to open the file.
    if ext in [".xls",".xlsx"]: #for Excel
        qc = pd.read_excel(CSV_converter)
    elif ext == '.csv': #for CSV
        qc = pd.read_csv(CSV_converter)
    elif ext == '.tsv': #for TSV
        qc = pd.read_csv(CSV_converter, sep="\t")
    
    return _convert_data(qc, channel)

def _convert_data(qc, channel):
    qc = qc.dropna() #removing nan, none from the data
    H = (qc['1H']*scale//12).astype(int) #scaling for proton, 0~12 ppm
    C = (qc['13C']*scale//240).astype(int) #scaling for carbon, 0~240 ppm
    if channel == 2: # When user chooses 'edited hsqc data' as input. the edited hsqc data should have three colums. '13C', '1H', and 'Intensitiy'
        T = qc['Intensity'] # To call intensity data
        #matrix for input
        mat = np.zeros((scale,scale,2), float)
        for j in range(len(qc)): 
            a, b = H.iloc[j].astype(int), C.iloc[j].astype(int)
            t = T.iloc[j]
            if 0 <= a < scale and 0 <= b < scale and t > 0:# + phase
                mat[b, scale-a-1,0] = 1 
            elif 0 <= a < scale and 0 <= b < scale and t < 0:# - phase
                mat[b, scale-a-1,1] = 1
    elif channel == 1: #When user chooses 'normal HSQC data' as input. 'Intensity' column isn't essential.
        mat = np.zeros((scale,scale,1), float)
        for j in range(len(qc)): 
            a, b = H.iloc[j].astype(int), C.iloc[j].astype(int)
            if 0 <= a < scale and 0 <= b < scale:
                mat[b, scale-a-1,0] = 1
    mat[:,0,:] = 1
    mat[:,127,:] = 1
    mat[0,:,:] = 1
    mat[127,:,:] = 1            
    return mat #it should be 128x128x1 for normal HSQC data, and 128x128x2 for edited HSQC data


#prediction function with two models
def predict_nmr(input_nmr_filename, channel, model): #filenamae, channel(normal/edited HSQC, 1/2)
    #To predict chemical fingerprints, molecular weights, compound class with probability, and glycoside checking with probability
    
    mat = CSV_converter(input_nmr_filename, channel) #output shape should be (128,128,1) or (128,128,2)
    ### Model Prediction
    ## Fingerprint and molecular weight prediction
    fingerprint_prediction, pred_MW, pred_class, pred_gly  = model.predict(mat.reshape(1,scale,scale,channel))

    return fingerprint_prediction, pred_MW, pred_class, pred_gly

# This takes the output of prediction from the keras model and post processes it
def refine_predict_nmr(fingerprint_prediction, pred_MW, pred_class, pred_gly):
    fingerprint_prediction = fingerprint_prediction[0].round() 
    fingerprint_prediction = np.where(fingerprint_prediction == 1)[0] #list of fingerprint
    pred_MW = pred_MW[0][0] # molecular weight
    
    
    ##Glycoside prediction
    pred_gly = pred_gly[0][1]
    ##Class prediction
    pred_class = pred_class[0]
    pred_class_index= sorted(range(len(pred_class)),key= lambda i: pred_class[i])[-3:] #top 5
    pred_class_index.reverse()
    pred_class_prob = pred_class[pred_class_index]

    #Return prediction results. Fingerprints, molecular weights, Compound class index, Compound class probability, and Glycoside probability
    return fingerprint_prediction, pred_MW, pred_class_index, pred_class_prob, pred_gly #array, array, array

def search_database(fingerprint_prediction, pred_MW, DB, mw=None, top_candidates=20):
    # To compare chemical fingerprints and molecular weights to find the molecules having similar properties (FP, MW)
    # Database structure, 2 must be the predictions for the DB, 3, must be the mass
    
    DB_len = len(DB)
    topK = np.full((DB_len,4), np.nan, dtype=object)
    for j in range(DB_len):
        if mw == None: #If we don't provide a user input mw, we should use the predicted MW
            DB_mw = DB[j][3]
            if abs(DB_mw-pred_MW)/(DB_mw) < 0.2 : #threshold is +- 20%
                DB_fp = DB[j][1]
                score = cosine(fingerprint_prediction,DB_fp)
                if score > 0.6:
                    topK[j] = DB[j][0], DB[j][2], score, DB_mw #Compound Name (list), SMILES (str), score, molecular weight of DB compounds 
        
        else: # If we provide a user input mw
            real_mw = DB[j][3]
            if abs(real_mw-pred_mw) <20:
                DB_fp = DB[j][1]
                score = cosine(fingerprint_prediction,DB_fp)
                if score > 0.6:
                    topK[j] = DB[j][0], DB[j][2], score, DB_mw #Compound Name (list), SMILES (str), score, molecular weight of DB compounds 
            

    
    #Saving the DB Search
    topK = pd.DataFrame(topK, columns = ['Name','SMILES','Cosine score','MW'] )
    topK = topK.dropna(how='all')
    topK = topK.sort_values(['Cosine score'], ascending = False)
    topK = topK.drop_duplicates(['SMILES'])
    topK = topK[:top_candidates] #remaining topN candidates
    topK = topK.fillna('No_name') #Time

    return topK


def _draw_search_image(output_nmr_image, mat, pred_class_index, pred_class_prob, pred_gly):
    mat[:,0,:] = 0
    mat[:,127,:] = 0
    mat[0,:,:] = 0
    mat[127,:,:] = 0  
    
    height, width = scale, scale
    plt.figure()
    ax = plt.axes()
    try:
        plt.imshow(mat[:,:,1], mpl.colors.ListedColormap([(0.2, 0.4, 0.6, 0),'blue']))
        plt.imshow(mat[:,:,0], mpl.colors.ListedColormap([(0.2, 0.4, 0.6, 0),'red']))
    except:
        plt.imshow(mat[:,:,0], mpl.colors.ListedColormap([(0.2, 0.4, 0.6, 0),'red']))
    ax.set_ylim(scale-1,0)
    ax.set_xlim(0,scale-1)
    print(pred_class_index, pred_class_prob,pred_gly)
    plt.axis()
    plt.xticks(np.arange(scale+1,step=scale/12), (list(range(12,0-1,-1))))
    plt.yticks(np.arange(scale+1,step=scale/12), (list(range(0,240+1,20))))
    ax.set_xlabel('1H [ppm]')
    ax.set_ylabel('13C [ppm]')
    plt.grid(True,linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    for n, idx in enumerate(pred_class_index):
        class_name = [i for i in index_super if index_super[i]==idx][0]
        plt.text(4, 8+8*n, f'{class_name} {round(float((100*pred_class_prob[n])),2)}%', ha='left', wrap=True, alpha=0.8)
    plt.text(4, 8+8*3, f'Glycoside {round(float((100*pred_gly)),2)}%', ha='left', wrap=True, alpha=0.8)
    #plt.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, hspace = 0, wspace = 0)
    plt.savefig(output_nmr_image, dpi=600)
    plt.close()


def search_CSV(input_nmr_filename, DB, model, channel, output_table, output_nmr_image, output_pred_fingerprint, mw=None, top_candidates=20): # i = CSV file name
    mat = CSV_converter(input_nmr_filename,channel)
    # Searching the molecules and drawing the HSQC spectra image with classificaiton results.
    # plotting and saving constructed HSQC images with predicted class
    ## image without padding and margin

    # This is the only time we need to call the model
    fingerprint_prediction, pred_MW, pred_class, pred_gly = predict_nmr(input_nmr_filename, channel, model)

    fingerprint_prediction, pred_MW, pred_class_index, pred_class_prob, pred_gly = refine_predict_nmr(fingerprint_prediction, pred_MW, pred_class, pred_gly)

    topK = search_database(fingerprint_prediction, pred_MW, DB, mw=mw, top_candidates=top_candidates)

    #Saving TopK
    topK.to_csv(output_table, index = None)
    
    #Saving the predicted Fingerprint
    open(output_pred_fingerprint, "w").write(json.dumps(fingerprint_prediction.tolist()))
    
    #Drawing constructed HSQC spectra
    #remove outline
    _draw_search_image(output_nmr_image, mat, pred_class_index, pred_class_prob, pred_gly)


def main():
    DB, index_super = load_db()
    model_1ch, model_2ch = load_models()

    parser = argparse.ArgumentParser(description='SMART Embedding')
    parser.add_argument('input_csv', help='input_csv')
    parser.add_argument('output_table', help='output_table')
    parser.add_argument('output_nmr_image', help='output_nmr_image')
    parser.add_argument('output_pred_fingerprint', help='output_pred_fingerprint')
    parser.add_argument('--molecular_weight', default=None, type=float, help='molecular_weight')
    parser.add_argument('--channel', default=1, type=int, help='HSQC type')
    args = parser.parse_args()
    
    if args.channel == 1:
        search_CSV(args.input_csv, DB, model_1ch, channel = args.channel, args.output_table, args.output_nmr_image, args.output_pred_fingerprint, mw=args.molecular_weight)
    elif args.channel == 2:
        search_CSV(args.input_csv, DB, model_2ch, channel = args.channel, args.output_table, args.output_nmr_image, args.output_pred_fingerprint, mw=args.molecular_weight)



if __name__ == "__main__":
    main()
