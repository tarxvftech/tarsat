run:
	python process.py nasabare.txt
update:
	#curl -4 https://live.ariss.org/iss.txt > iss.txt
	#python process.py iss.txt
	#https://www.amsat.org/tle/current/nasabare.txt
	curl -4 https://www.amsat.org/tle/current/nasabare.txt > nasabare.txt
clean:
	-rm -r generated
	mkdir generated
