# pgen-tmle-genetic-interactions
Run Interaction Analyses from PGEN files

## Building the Docker Image

```bash
docker build --platform linux/amd64 -t olivierlabayle/pgen-tmle-interactions:main -f docker/Dockerfile .
```

## Running the Docker Image


```bash
docker run -it --platform linux/amd64 olivierlabayle/pgen-tmle-interactions:main /bin/bash
```