FROM haoencui/baker

RUN mkdir /home/ica

RUN R -e "options(repos = \
   list(CRAN = 'http://cran.ma.imperial.ac.uk')); \
   install.packages('R2jags'); \
   install.packages('BCEA')"

