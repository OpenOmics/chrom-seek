## Steps for Building Docker Images

Directly below are instructions for building a base image for the `cfchip_toolkit` using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=cfchip_toolkit:v0.2.0 .

# Testing, take a peek inside
docker run -ti cfchip_toolkit:v0.2.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag cfchip_toolkit:v0.2.0 skchronicles/cfchip_toolkit:v0.2.0
docker tag cfchip_toolkit:v0.2.0 skchronicles/cfchip_toolkit         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/cfchip_toolkit:v0.2.0
docker push skchronicles/cfchip_toolkit:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan cfchip_toolkit:v0.2.0
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
