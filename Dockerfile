FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.5

WORKDIR /build
COPY src .
RUN make

WORKDIR /app
ENV PATH=$PATH:/app
RUN cp /build/covcollect /app
