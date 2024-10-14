# antarcticaRadioNeutrino

This is a simple set of examples showing how we can read information stored in HDF5 format describing the kind of data that might be acquired on a radio balloon experiment searching for neutrinos in Antarctica

## Ballon Plotters

There are three simple notebooks included here:
    1. ballonHdf5HeaderReader.ipynb  -- Which shows how we can read the small amount of information included in the event header (along with position and attitude from the GPS)
    2. ballonMapPlotter.ipynb -- Which shows how we can plot some of the positin variables on an Antarttica projection
    3. ballonHdf5EventReader.ipynb -- Shows how to read the large amount of waveform data in a run and construct simple images.

