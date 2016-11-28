# ![fourier-being](logo.png)
<img src="https://ci.appveyor.com/api/projects/status/tjvo1hcohnrkx58g/branch/master?svg=true" align="center" />

# Least Squared Error Based FIR Filters
> A very flexible graphical interface to design Least Square Error Based FIR Filters and generate fixed point coefficients 

## Getting Started
### Install Python on Ubuntu
```
http://askubuntu.com/questions/101591/how-do-i-install-python-2-7-2-on-ubuntu
```
#### NOTE: Follow similar instructions on https://www.python.org/ to get python on windows or macOS
### Get the repositiory
```
username@ubuntu:~$ git clone https://github.com/fourier-being/Least-Squared-Error-Based-FIR-Filters.git
```
### Go to the Repo Directory
```
username@ubuntu:~$ cd Least-Squared-Error-Based-FIR-Filters/
```
### Run the Python script
```
username@ubuntu:~/Least-Squared-Error-Based-FIR-Filters$ python LSDesignAdvanced.py
```
> After running the python script, you should see the following GUI :

<img src="figure_1.png" align="center" />

#### NOTE: There are sliding bars at the bottom of the plot window corresponding to different parameters. You can adjust the sliding bars according to your requirement and then click on "Generate" button to generate the coefficients.

## Description
This tool currently implements only Low-Pass FIR Filter. I am working to add other filters - Bandpass, High Pass and Multiple Band Pass.
Future scope includes a python package which would have all these FIR filters APIs and Coefficients generated using this GUI based tool. Least-Squared Error based FIR Filters basically minimizes the squared error between desired filter response and amplitude response of the FIR filter. I am taking Type-I filter here since it is Low-Pass.

## License

MIT © [fourier-being]
