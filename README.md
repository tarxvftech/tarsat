

## TODO
Define and implement useful APIs for pre-calculating passes, storing them, and displaying them.
Offer a semi-standard API for predicting passes.
Draw doppler curves and prove they are accurate.

Make test cases of specific times and ephemeris to allow validating different algorithms.
Characterize ephemerides over time and make a blog post investigating them.
    Since the ephemerides get updated, and should be accurate over time
    if not precise, we can show the differences over time and show how
    each N-M ephemeride compares to the ones after it.
    One of the links below has some discussion of this, too.

For a simple test
    load a particular TLE for a satellite from file
        make this easy!
    pick a specific startjd
    output ra/dec/...alt? whatever it is
    for a location, plot the pass
    for that location, plot the doppler curve for a theoretical transmitter



Links
https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=2075&context=smallsat
https://conference.sdo.esoc.esa.int/proceedings/neosst1/paper/435/NEOSST1-paper435.pdf




## Ideas:
* framebuffer / sixel drawing
    https://github.com/hackerb9/lsix
    https://www.nongnu.org/fbi-improved/
