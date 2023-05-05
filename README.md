# FAST-CASPI

To run need to download fftw: https://www.fftw.org/download.html

Need to configure fftw by using command in terminal "./configure --enable-threads", "make", and "make install" in order.

Then run CASPI on terminal using makefile.

PLEASE NOTE: For optimal time FFTW needs to ran first to write fftw plans to fftw wisdom, then ran again to read those plans. This is currently commented out, so run time is slower if ran. 
