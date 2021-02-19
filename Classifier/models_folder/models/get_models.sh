rm CHANNEL1 CHANNEL2-r
rm *hdf5

wget --output-document=models.zip "https://www.dropbox.com/sh/ndtdcs212h7gvwn/AACswt63ScOxrWN6auK3V-Cca?dl=1"
unzip models.zip
mkdir -p CHANNEL1/000001
mkdir -p CHANNEL2/000001

python ./convert_keras_to_tf.py "(020321)SMART3_v3_1ch_multitask_final.hdf5" CHANNEL1/000001
python ./convert_keras_to_tf.py "(020321)SMART3_v3_2ch_multitask_final.hdf5" CHANNEL2/000001

