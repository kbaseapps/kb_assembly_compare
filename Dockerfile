FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install coverage

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade


#### Install MUMmer4
##
RUN mkdir -p /kb/module
WORKDIR /kb/module
# YAGGO
RUN curl -s https://codeload.github.com/gmarcais/yaggo/tar.gz/v1.5.10 > yaggo.v1.5.10.tar.gz && \
    tar xfz yaggo.v1.5.10.tar.gz && \
    cd yaggo-1.5.10 && \
    make DEST=/usr/local/bin && \
    cp ./yaggo /usr/local/bin
# MUMmer4
#RUN git clone https://github.com/mummer4/mummer && \
#    cd mummer && \
##    mkdir aux_bin && \
##    make check && \
#    make install CPPFLAGS="-O3 -DSIXTYFOURBITS"
##    mummerpath=`pwd` && \
##    cd .. && \
##    export PATH=$mergerpath:$mummerpath:$PATH

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
