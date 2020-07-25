# toolkit
All these are wheels that I can hardly remember how many times that I reiteratively write them.

xrd_main.py
-----------
A primary code that simulate x ray diffraction, without nearly any acceleration modification on code, except smooth function built in it.
To be honest in the past I just run a smooth with time complexity of O(N**4) which seems a bit silly now but, I used to have much time to wait and stare at the screen...

features:
1. totally arbitrary wavelength option
2. totally arbitrary diffraction angle setting. but use with care, I think time complexity will be scaring.
3. only single-element structures are supported in present version. Difference between kinds of atoms will reflect in diffraction coefficients but it seems that I forget where I stored that paper parameterizing coefficients. :( 
   I will find it.
4. angle resoluted diffraction coefficients is implemented, but coefficients need to be supplied manually in this version. I will improve it soon, also.
5. supercell is not implemented yet, but I have completed this function in other seperated code. this xrd simulation code will be deprecated if I finish combining all functions together.
6. Simple cutoff kind method denoise is implemented in this code, but I will find ways to modify it and maybe conserve more data points if necessary.
7. simple gaussian smooth is implemented. I also get inspired when writing this, i.e., use cutoff when smoothing the whole 2d picture.

coordsRead.py
-------------
it will read standard-formatted *.xyz files, convert coordinates into list, there is also option that can convert str-type data into float, which will be beneficial for latter processing.

supercell.py
------------
duplicate a cell any times you want. six cell parameters are required, also output of coordsRead.py can be directly used in parameter list of this subroutine.
this subroutine will return a duplicated, formated list, containing 3 or 4 columns depending on data input.

smooth.py
---------
1d and 2d gaussian function with or without cutoff acceleration method is implemented. it is suitable for ALL data that does not exceeds two dimensions, which is aim to improve
visual experience at the beginning.
Fermi-dirac function can not work normally till now, I will check.

%Recent plan
I will implement FFT to accelerate diffraction module, and maybe in some day I will integrate all them together, yielding a home-made simulation package.
:) Haha...
