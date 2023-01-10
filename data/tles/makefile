run:
	python process.py nasabare.txt
update:
	#curl -4 https://live.ariss.org/iss.txt > iss.txt
	#python process.py iss.txt
	#https://www.amsat.org/tle/current/nasabare.txt
	curl -4 https://www.amsat.org/tle/current/nasabare.txt > nasabare.txt
	./rename.sh nasabare.txt
	#need historical data. Can start dumping nasabare to disk/git.
	#can also fetch from
	#https://celestrak.org/NORAD/archives/
	#https://www.space-track.org/auth/login
	#https://www.amsat.org/tle/
clean:
	-rm -r generated
	mkdir generated
