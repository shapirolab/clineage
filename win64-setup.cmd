@echo off

echo ===============================
echo Installing lxml
easy_install lxml==2.3

echo ===============================
echo Installing MySQL-python
easy_install MySQL-python==1.2.5

echo ===============================
echo Installing NumPy
easy_install http://pkgs.10x.org.il/amd64/numpy-MKL-1.8.1.win-amd64-py2.7.exe

echo ===============================
echo Installing biopython
easy_install http://pkgs.10x.org.il/amd64/biopython-1.63.win-amd64-py2.7.exe

echo ===============================
echo Installing pycrypto
easy_install http://www.voidspace.org.uk/downloads/pycrypto26/pycrypto-2.6.win-amd64-py2.7.exe

echo ===============================
echo Installing paramiko
pip install paramiko==1.10.2