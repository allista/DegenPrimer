#!/bin/bash

inst()
{
    package=$(basename $(pwd))
    echo Removing /usr/lib/python2.7/site-packages/$package*
    sudo trash-put /usr/lib/python2.7/site-packages/$package*
    echo Removing /usr/local/lib/python2.7/site-packages/$package*
    sudo trash-put /usr/local/lib/python2.7/site-packages/$package*
    python setup.py build && sudo python setup.py install
}


pushd BioUtils
inst

popd
inst

pushd ../DegenPrimerGUI
inst
