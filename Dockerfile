FROM lucas501/cytonorm_docker:latest
#FROM tercen/runtime-r40-slim:4.0.4-0

RUN git clone https://github.com/tercen/cytonorm_operator.git

WORKDIR /operator/cytonorm_operator

RUN echo 0.0.11 && git pull
RUN git checkout 0.0.11

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT [ "R","--no-save","--no-restore","--no-environ","--slave","-f","main.R", "--args"]
CMD [ "--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
