# 1. Not Another Variant Annotator

Another variant annotator that leverages the ExAC server to 
annotate ones vcf. 

*Caveats*

This has only been tested on FreeBayes output for two triploid samples.

It only works on the human genome assembley GR37 as 
that is all ExAC support. Also their version of VEP, which performs the annotations for the impact of the variant is rather outdated.


## 1.1. Getting Started

There isn't much required in terms of pre-requistes, besides an internet connection. The script does support the use of snpEff, however one must already have a local copy of the software installed. This can be found here https://pcingola.github.io/SnpEff

## 1.2. Install

clone the repo to your local machine
cd into the directory of the cloned repo and install with:
`python setup.py install`

## 1.3. Running the tests

Limited unit testing has been made for this code using a combination of
the hypothesis library to perform property-based testing and Contracts.

Example command to run the test:

`pytest /path_to_repo/src/tests`


## 1.4. Built With

* [py3.8](https://www.python.org/downloads/release/python-380/) - python
* [pandas](https://pandas.pydata.org/) - python library for data manipulation
* [hypothesis](https://hypothesis.works/) - python library for generating test data
* [dpcontracts](https://github.com/deadpixi/contracts) - python library for setting contracts

## 1.5. Versioning

No verision software is used at this time

## 1.6. How to use

A simple tool to annotate one's vcf.

Example command to run the annotator script with all annotation functions:

`NAVANN --vcf ~/example/Challenge_data.vcf --verbose --outpath /path/to/outfolder/Challenge_data.ann.vcf --all`

Run the annotator script with only adding vaf to the format field:

`NAVANN --vcf ~/example/Challenge_data.vcf --verbose --outpath /path/to/outfolder/Challenge_data.ann.vcf --vaf`


## 1.7. Authors

* **Andrei Rajkovic** - *Initial work* 

## 1.8. License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details