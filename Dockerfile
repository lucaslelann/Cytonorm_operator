#FROM tercen/runtime-r40:4.0.4-1
FROM tercen/runtime-r40-slim:4.0.4-0


USER root
WORKDIR /operator/

RUN apt-get update
RUN apt-get install -y r-cran-tcltk2

RUN git clone https://github.com/tercen/cytonorm_operator.git

WORKDIR /operator/cytonorm_operator 

RUN git checkout master
RUN echo 0.0.2 && git pull
RUN git checkout 0.0.2

RUN R -e "renv::consent(provided=TRUE);renv::restore(confirm=FALSE)"

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT [ "R","--no-save","--no-restore","--no-environ","--slave","-f","main.R", "--args"]
CMD [ "--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
