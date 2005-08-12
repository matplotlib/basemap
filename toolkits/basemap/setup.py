from distutils.core import setup, Extension
import sys

deps = ['src/PJ_aeqd.c','src/PJ_gnom.c','src/PJ_laea.c', \
	'src/PJ_mod_ster.c','src/PJ_nsper.c','src/PJ_nzmg.c','src/PJ_ortho.c','src/PJ_stere.c', \
	'src/PJ_sterea.c','src/PJ_aea.c','src/PJ_bipc.c','src/PJ_bonne.c','src/PJ_eqdc.c',\
	'src/PJ_imw_p.c','src/PJ_krovak.c','src/PJ_lcc.c','src/PJ_mpoly.c','src/PJ_poly.c',\
	'src/PJ_rpoly.c','src/PJ_sconics.c','src/PJ_cass.c','src/PJ_cc.c','src/PJ_cea.c',\
	'src/PJ_eqc.c','src/PJ_gall.c','src/PJ_labrd.c','src/PJ_lsat.c','src/PJ_merc.c',\
	'src/PJ_mill.c','src/PJ_ocea.c','src/PJ_omerc.c','src/PJ_somerc.c','src/PJ_tcc.c',\
	'src/PJ_tcea.c','src/PJ_tmerc.c','src/PJ_airy.c','src/PJ_aitoff.c','src/PJ_august.c',\
	'src/PJ_bacon.c','src/PJ_chamb.c','src/PJ_hammer.c','src/PJ_lagrng.c','src/PJ_larr.c',\
	'src/PJ_lask.c','src/PJ_nocol.c','src/PJ_ob_tran.c','src/PJ_oea.c','src/PJ_tpeqd.c',\
	'src/PJ_vandg.c','src/PJ_vandg2.c','src/PJ_vandg4.c','src/PJ_wag7.c','src/PJ_lcca.c',\
	'src/PJ_geos.c','src/PJ_boggs.c','src/PJ_collg.c','src/PJ_crast.c','src/PJ_denoy.c',\
	'src/PJ_eck1.c','src/PJ_eck2.c','src/PJ_eck3.c','src/PJ_eck4.c','src/PJ_eck5.c',\
	'src/PJ_fahey.c','src/PJ_fouc_s.c','src/PJ_gins8.c','src/PJ_gn_sinu.c','src/PJ_goode.c',\
	'src/PJ_hatano.c','src/PJ_loxim.c','src/PJ_mbt_fps.c','src/PJ_mbtfpp.c',\
	'src/PJ_mbtfpq.c','src/PJ_moll.c','src/PJ_nell.c','src/PJ_nell_h.c','src/PJ_putp2.c',\
	'src/PJ_putp3.c','src/PJ_putp4p.c','src/PJ_putp5.c','src/PJ_putp6.c','src/PJ_robin.c',\
	'src/PJ_sts.c','src/PJ_urm5.c','src/PJ_urmfps.c','src/PJ_wag2.c','src/PJ_wag3.c',\
	'src/PJ_wink1.c','src/PJ_wink2.c','src/PJ_latlong.c','src/PJ_geocent.c',\
	'src/aasincos.c','src/adjlon.c','src/bch2bps.c','src/bchgen.c','src/biveval.c',\
	'src/dmstor.c','src/mk_cheby.c','src/PJ_auth.c','src/PJ_deriv.c','src/PJ_ell_set.c',\
	'src/pj_ellps.c','src/PJ_errno.c','src/PJ_factors.c','src/PJ_fwd.c','src/PJ_init.c',\
	'src/pj_inv.c','src/PJ_list.c','src/PJ_malloc.c','src/PJ_mlfn.c','src/PJ_msfn.c',\
	'src/pj_open_lib.c','src/PJ_param.c','src/PJ_phi2.c','src/PJ_pr_list.c','src/PJ_qsfn.c',\
	'src/pj_strerrno.c','src/PJ_tsfn.c','src/PJ_units.c','src/PJ_zpoly1.c','src/rtodms.c',\
	'src/vector1.c','src/PJ_release.c','src/PJ_gauss.c','src/nad_cvt.c','src/nad_init.c',\
	'src/nad_intr.c','src/emess.c','src/PJ_apply_gridshift.c','src/PJ_datums.c',\
	'src/pj_datum_set.c','src/PJ_transform.c','src/geocent.c','src/PJ_utils.c',\
	'src/pj_gridinfo.c','src/PJ_gridlist.c','src/proj4.c']


extensions = [Extension("proj4",
                        deps,
                        include_dirs = ['src'],)]

datadir ='share/basemap-py'+repr(sys.version_info[0])+repr(sys.version_info[1])

setup(
  name              = "basemap",
  version           = "0.4",
  description       = "Plot data on map projections with matplotlib",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  data_files        = [(datadir,['data/countries_c.txt','data/states_c.txt','data/countries_l.txt','data/states_l.txt','data/gshhs_c.txt','data/gshhs_l.txt','data/countries_i.txt','data/states_i.txt','data/gshhs_i.txt'])],
  packages          = ['matplotlib/toolkits','matplotlib/toolkits/basemap'],
  package_dir       = {'':'lib'},
  ext_modules       = extensions)
