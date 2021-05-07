## SMART 3

### Checking Model Metadata

We pass through tensorflow serving at this url:

```/model/metadata```

If the model input names change, then we need to change it in the code

### APIs

Classify programmatically 

```/api/smart3/search```

You can put in your peaks as a json list of dicts, with 1H,13C as headers. 

## License

The license as included for the software is MIT. Additionally, all data, models, and ontology are licensed as [CC0](https://creativecommons.org/share-your-work/public-domain/cc0/).

## Privacy

We try our best to balance privacy and understand how users are using our tool. As such, we keep in our logs which structures were classified but not which users queried the structure. 
