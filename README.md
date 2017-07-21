# UVic Earth System Climate Model Oregon State University CEOAS Andreas Schmittner Group 

Workspace for modified code of [University of Victoria Earth Model](http://climate.uvic.ca/model/) by Dr. Andreas Schmittner's
physical oceanography group at the College of Environment, Ocean, and Atmospheric Science at Oregon State University.


## Getting Started

The master branch is the 'original' code from the University of Victoria model. You'll need to check out the code and compile. Branches correspond modifications of code to align with different research goals.

### Prerequisites

```
Give examples
```

### Installing

* First, update `mobi1.9.q` so that the contents reflect your current working directory. This is the file that you submit to the job scheduler.

```
#!/bin/csh
# BE SURE To UPDATE YOUR WORKING DIRECTORY BELOW
#$ -e /home/changeme/UVic2.9/MOBI1.9
#$ -o /home/changeme/UVic2.9/MOBI1.9
cd /home/changeme/UVic2.9/MOBI1.9
# the next line is only needed on the magellan linux cluster at CEOAS
source /share/apps/lf6481/csh_laheyfort_setup
time ./UVic_ESCM > pr
```

Assuming everything compiled correctly, you should be able to submit a job:

```
$ qsub mobi1.9.q 
```

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **University of Victoria** - *Initial work* - [UVic](http://climate.uvic.ca/model/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc


