# PhyloPro

PhyloPro is a phylogenetic profile pipeline that takes advantage of the MiST and SeqDepot database.

### Get Started

#### Step 1 - Kickstart the software

``` bash
PhyloPro.py --init ProjectName
```

make sure to not pick a ProjectName with special characters as they will become part of the filename used by **PhyloPro**.

This will create the system of directories and a config file `phylopro.ProjectName.cfg.json` at the main project directory.

#### Step 2 - Edit the config file

Open the `phylopro.ProjectName.cfg.json` in a text editor and edit.

#### Step 3 - Continue PhyloPro

``` bash
PhyloPro.py --continue ProjectName
```




### The pipeline
The pipeline is serial but you can restart the task from any point in case of needed intervention by passing the flag. If for same reason the pipeline fails, it will restart from where stopped which is recorded in the file ```restart.phylopro.json```. Below is the list of steps in the pipeline and the flags to be passed in case to overwrite the restart.

#### Pipeline

Prepare files and databases to be used in the future

| Flag name           | Default Argument value | Description |
|---------------------|--------------| ---------------|
| --init 			  | ProjectName | Build the local directory system and writes the local configuration file: phylopro.local.cfg.json|
| --fetchProtFams     | ProjectName | Acquire and parse info about protein families from SeqDepot |
| --fetchGenInfo      | ProjectName | Acquire and process genome information from MiST3 |
| --mkFastaFiles      | Make relevant fasta files |
| --filterByGen       | Filter fasta files to only contain sequences from the selected genomes |
| --trimSeqs          | Trim sequences based on HMM models |

### Requisites
To use PhyloPro you need a few items:

1. The MiST3 API in python	

### Configuration
PhyloPro requires two config files. A global one at the root of the PhyloPro directory ```phylopro.cgf.json``` and a local one in the root directory of the job ```phylopro.local.cfg.json```

#### The configuration files
The configuration files are JSON formatted and contain information necessary to run the pipeline. 

##### phylopro.cfg.json

``` javascript
[ { mistdir : "" // Directory to where all mist files are}]
```

##### phylopro.local.cfg.json

``` json
[ { "genomes" : 
		[ 000, 001, 002 , ... , N ], // List of mistID of genomes for analysis
	
	"ProtFamDef" : { // Define protein family using results by tools in
		"SeqDepot"	:
		[ 	{ "name" : "Component1",
			  "group" : "Group name 1"
			  "pfam29" : { 	"in" : ["PfamDomain1, PfamDomain2", ... ], 
						 	"out": ["PfamDomain3, PfamDomain4", ... ]
				}
			},
			{	"name" : "Component2"
				"pfam29" : {
					"in"  : ["PfamDomain1, PfamDomain2", ... ], 
					"out" : ["PfamDomain3, PfamDomain4", ... ]
				}
			},
		],
		"CustomHMM" : 
		[	{ "name" : "Component1",
			  "HMMfile" : "path_to_file_HMM1"
			},
			{ "name" : "Component2",
			  "HMMfile" : "path_to_file_HMM2"
			}
		
	]
```


### Directory scheme

#### source

./MiSTFiles/	--			where the MiST fasta files will exists

#### local



### More about the pipeline

#### Select the Protein family definition

The selection of protein family will be done now based on the architecture depending on the tools ran by SeqDepot.


#### Example of Local config file:

``` json
{
  "ProtFamDef": {
    "SeqDepot": [
      {
        "group": "Adaptors", 
        "name": "CheW", 
        "pfam29": {
          "out": [
            "HATPase_c", 
            "Response_reg"
          ], 
          "in": [
            "CheW"
          ]
        }
      }, 
      {
        "group": "Adaptors", 
        "name": "MCPs", 
        "pfam29": {
          "in": [
            "MCPsignal"
          ]
        }
      }, 
      {
        "group": "Histidine Kinase", 
        "name": "CheA", 
        "agfam1": {
          "in": [
            "HK_CA:Che"
          ]
        }
      }, 
      {
        "group": "Histidine Kinase", 
        "name": "CheA", 
        "pfam29": {
          "in": [
            "CheW", 
            "HATPase_c"
          ]
        }
      }
    ]
  }, 
  "ProjectName": "Vibrio", 
  "stage": "init", 
  "genomes": [
    479, 
    319
  ], 
  "history": {
    "mkFastaFiles": "2016-05-06T20:32:50.324466", 
    "fetchProtFams": "2016-05-06T20:32:49.601795", 
    "init": "someday", 
    "fetchGenInfo": "2016-05-06T20:32:50.322592"
  }
}```
