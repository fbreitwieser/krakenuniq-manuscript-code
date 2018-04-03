
all: install

krakenhll:
	git clone https://github.com/fbreitwieser/krakenhll

install:
	krakenhll/install_krakenhll.sh `pwd`/install