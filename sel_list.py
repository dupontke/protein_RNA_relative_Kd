##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

sel = []
sel.append(['resid 50 and not backbone','resid 28 and not backbone','R50E28'])
sel.append(['resid 50 and not backbone','resid 104 and not backbone','R50E104'])
sel.append(['resid 18 and not backbone','resname GTP and name N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4','R18gtp_base'])
sel.append(['resid 50 and not backbone','resname GTP and name O1G PG O2G O3G','R50gtp_gamma_phosphate'])
