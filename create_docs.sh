sudo pip uninstall sibreg
sudo python setup.py install
sphinx-apidoc -f -o docs sibreg
make -C docs clean html 