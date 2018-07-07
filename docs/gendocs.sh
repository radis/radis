rm -rf source/ 
sphinx-apidoc -f -o source/ ../radis
make clean
make html
