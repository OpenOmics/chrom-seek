## Steps for Building Docker Images

Directly below are instructions for building a base image for the `chrom_seek_python` using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=chrom_seek_python:v0.1.0 .

# Testing, take a peek inside
docker run -ti chrom_seek_python:v0.1.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag chrom_seek_python:v0.1.0 org_account/chrom_seek_python:v0.1.0
docker tag chrom_seek_python:v0.1.0 org_account/chrom_seek_python         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push org_account/chrom_seek_python:v0.1.0
docker push org_account/chrom_seek_python:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan chrom_seek_python:v0.1.0
```

> **Please Note**: Any references to `org_account` should be replaced your username if you would also like to push the image to a non-org account.