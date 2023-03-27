IMAGE_NAME=svpop
VERSION=0.1

TAG1=wharvey31/$(IMAGE_NAME):$(VERSION)
TAG2=wharvey31/$(IMAGE_NAME):latest

all: build push

build:
	sudo docker build -t $(TAG1) -t $(TAG2) .

push:
	sudo docker push $(TAG1)
	sudo docker push $(TAG2)
