
all: install

krakenuniq:
	git clone https://github.com/fbreitwieser/krakenuniq

install:
	krakenuniq/install_krakenuniq.sh `pwd`/install