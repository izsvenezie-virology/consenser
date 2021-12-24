# Consenser
Creates a consensus sequence from a vcf file and a reference.

## Installation:
To install consenser:
```
cd /path/where/download
git clone https://gitlab.izsvenezie.it/EdoardoGiussani/consenser
cd consenser
python -m pip install .
```
To test if it is installed correctly type 
```
consenser --help
```
and the help page should appear on your screen. 
If this not happens you can write me threatening letters.


### Installation issues
If you have some problems with SSL Certificates when trying to clone the git repository type:
```
git config --global http.sslVerify false
```
Don't forget to reactivate SSL Certificates once you cloned it or hackers will steel all your passwords:
```
git config --global http.sslVerify true
```

## Usage
To use consenser type:
```
consenser [OPTIONS] 
```
| Option                            | Description                                   | Default           |
| --------------------------------- | --------------------------------------------- | ----------------- |
|```-c```/```--coverage FILENAME``` | The coverage file to use, a tab separated file with columns 'Chromosome', 'Position', 'Coverage'. To use the stdin instead of a file use '-'.                 | REQUIRED          |
|```-r```/```--reference FILENAME```| The reference file in Fasta format. To use the stdin instead of a file use '-'.                                                                                | REQUIRED          |
|```-v```/```--vcf FILENAME```      | The vcf file. Actually works with the Lofreq vcf format. To use the stdin instead of a file use '-'.                                                                     | REQUIRED          |
|```-o```/```--output FILENAME```   | The produced consensus file. To print the consensus to stdout use '-'.                                                                                | REQUIRED          |
|```-d```/```--deg FLOAT FLOAT```   | Upper and lower degeneration limits.          | -                 |
|```-n```/```--no-deg FLOAT```      | Minimum frequency to report a mutation.       | -                 |
|```-m```/```--min-cov INT```       | The minimum coverage to report a nucleotide.  | 10                |
|```-w```/```--width INT```         | Maximum width of the Fasta lines.             | 70                |
|```-s```/```--split PATH```        | Splits the consensus by chromosome.           | -                 |
|```-a```/```--alter-names TEXT```  | Change chromosomes name.                      | -                 |
|```--indels-lim FLOAT```           | Set minimum limit to consider an indel.       | 50.0              |

Options ```--deg``` or ```--no-deg``` are mutually exclusive and required.

## Degenerated consensus
If you want a consensus with IUPAC degenerated bases you shoud use the ```--deg``` option.
This option requires 2 float numbers. 
These indicates the lower and upper limit of nucleotides frequency (in percentage).
Nucleotides are used with the following rules:
* If the nucleotide frequence is lower than the lower limit, is discarded.
* If the nucleotide frequence is higher than the upper limit, is reported in the consensus.
* If the nucleotide frequence is between the two limits included, is used in the degeneration.

You shuold use a specular interval (eg. ```--deg 25 75```), you are free to do what you want, but I can't assure the results are correct.

In the degeneration conversion is also considered the reference base (if frequency is higher than lower limit).

## Non degenerated consensus
If you want a consensus with only A, C, T, G and N use the option ```--no-deg```. 
This option requires a limit (in percentage).

* If the nucleotide frequency is greater than the limit the nucleotide is reported in the consensus.
* If none of the nucleotides is greater than the limit the reference nucleotide is reported in the consensus.

The suggested value is ```--no-deg 50```. You can use an higher percentage, but not a lower percentage because it's nonsense and the results will be probably wrong.

## Minimum coverage
With the option ```--min-cov``` you can set the threshold for mask the low covered positions.
Positions with a coverage below the treshold is masked with the base "N". 

To skip this check use ```--min-cov 0```.


