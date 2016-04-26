# PhyloPro

PhyloPro is a phylogenetic profile pipeline that takes advantage of the MiST and SeqDepot database.

### Get Started

``` bash
PhyloPro.py --init
```

### The pipeline
The pipeline is serial but you can restart the task from any point in case of needed intervention by passing the flag. If for same reason the pipeline fails, it will restart from where stopped which is recorded in the file ```restart.phylopro.json```. Below is the list of steps in the pipeline and the flags to be passed in case to overwrite the restart.

#### Pre-pipeline

Prepare files and databases to be used in the future

| Flag name           | Default Argument value | Description |
|---------------------|--------------| ---------------|
| --init 			  | ./  | Build the local directory system and writes the local configuration file: phylopro.local.cfg.json|
| --build-MiST-fasta  | ./MiSTFiles/ | Required to use your own HMM to select sequences from protein families. It will automatically update the global config file.   |

#### Main Pipeline

| Flag name           | Description |
|---------------------| ---------------|
| --files			  | Make the relevant files |
| --man-files	      | Manipulate fasta files: Trim |


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

``` javascript
[ { 
	ProtFamDef :  // Define protein family using results by tools in
		SeqDepot	:
		[ 	{ name : 'Component1',
			  pfam29 : { in : ['PfamDomain1, PfamDomain2', ... ], 
						 out: ['PfamDomain3, PfamDomain4', ... ]
				}
			},
			{	name : 'Component2'
				pfam29 : {
					in : ['PfamDomain1, PfamDomain2', ... ], 
					out: ['PfamDomain3, PfamDomain4', ... ]
				}
			},
		],
		CustomHMM : 
		[	{ name : 'Component1',
			  HMMfile : 'path_to_file_HMM1'
			},
			{ name : 'Component2',
			  HMMfile : 'path_to_file_HMM2'
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