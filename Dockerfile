ARG platform=linux/amd64
FROM --platform=$platform ubuntu:24.04
LABEL authors="embarkvet"

WORKDIR /build

RUN apt update \
  && apt install -y \
    build-essential \
    g++ \
    gcc

COPY . /build
