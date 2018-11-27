FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update

# To install all the dependencies
RUN apt-get update && apt-get install -y build-essential wget make curl unzip python && \
    apt-get install -y r-base r-cran-gplots

# Install pandas
#RUN pip install pandas

# Install X spoof execution wrapper (for plotting)
RUN apt-get -y install xvfb

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
