import astropy.units as u
from astropy.coordinates import SkyCoord

def catalog_match(ra1, dec1, ra2, dec2):
    c1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    c2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
    idx, d2d, _ = c1.match_to_catalog_sky(c2)

    return idx, d2d

def to_ds9_region(ra, dec, text=None, filename="ds9.reg"):
    with open(filename,'w+',newline='') as txtfile:
        txtfile.write('# Region file format: DS9 version 4.1\n')
        txtfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        txtfile.write('fk5\n')
        if text is None:
            for _ra, _dec in zip(ra, dec):
                txtfile.write(f"point({_ra},{_dec}) # point=diamond width=2 color=cyan \n")
        else:
            for _ra, _dec, _text in zip(ra, dec, text):
                txtfile.write(f"point({_ra},{_dec}) # point=diamond width=2 color=cyan")
                txtfile.write(" text={" + str(_text) + "}\n")