'''
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
'''
import string, os, math

f = open("input1","r")   # input file for the python program
ff = open("input","w")   #input file for the moldft program

#get the value of the electric field from the input file

def get_field(fobj):
    f = 0.0
    lines = fobj.readlines()
    for line in lines:
        ind = string.find(line,"Field")
        if ind != -1:
            f = float(line.split()[1])
            return f
    fobj.close()
#get the energy of the molecule in electric field
def get_energy(field):
    if field ==" ":
        field = 'nf'
    else:
        field = field
    
    f = open("staticdir/static_%s.out"%field,"r")
    while 1:
        line = f.readline()
        ind = string.find(line,'Summary')
        if ind !=-1:
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            line = f.readline()
            line1 = line.split()
            return float(line1[1])
        if not line:
            break
#function to call moldft and compute dipole 
def process_file(field,List):
    if field ==" ":
        field = 'nf'
    else:
        field = field
    f = open("input1","r")   # input file for the python program
    ff = open("input","w")   #input file for the moldft program
    lines  = f.readlines()
    for line in lines:
        ind = string.find(line,"field_comp")
        if ind != -1:
            line = line.strip("field_comp\n")+ "%s"%field +"\n"
        ff.write(line)
    f.close()
    ff.close()
    #place molecule in field
    #os.mkdir('staticdir_1au')
    os.system("./moldft>staticdir/static_%s.out"%field)
    ffx = open("staticdir/static_%s.out"%field,"r")
    lins = ffx.readlines()
    for line in lins:
        ind1 = string.find(line, "mux")
        ind2 = string.find(line, "muy")
        ind3 = string.find(line, "muz")
        if ind1 !=-1:
            muxfx = float(line.split()[1])
            List.append(muxfx)
        if ind2 !=-1:
            muyfx = float(line.split()[1])
            List.append(muyfx)
        if ind3 !=-1:
            muzfx = float(line.split()[1])    
            List.append(muzfx)
            if not line:
                break
    ffx.close()

# molecule in fx field
field = "fx"
fxdip = []
fx_field = process_file(field,fxdip)
uxfx = float(fxdip[0])
uyfx = float(fxdip[1])
uzfx = float(fxdip[2])
Efx = get_energy(field)
# molecule in fy field
fyfield ="fy"
fydip =[]
fy_field = process_file(fyfield,fydip)
uxfy = float(fydip[0])
uyfy = float(fydip[1])
uzfy = float(fydip[2])
Efy = get_energy(fyfield)
# molecule in fz field
fzfield ="fz"
fzdip =[]
fz_field = process_file(fzfield,fzdip)
uxfz = float(fzdip[0])
uyfz = float(fzdip[1])
uzfz = float(fzdip[2])
Efz = get_energy(fzfield)
# molecule in -fx field
mfxfield ="mfx"
mfxdip =[]
mfx_field = process_file(mfxfield,mfxdip)
uxmfx = float(mfxdip[0])
uymfx = float(mfxdip[1])
uzmfx = float(mfxdip[2])
Emfx = get_energy(mfxfield)
# molecule in -fy field
mfyfield ="mfy"
mfydip =[]
mfy_field = process_file(mfyfield,mfydip)
uxmfy = float(mfydip[0])
uymfy = float(mfydip[1])
uzmfy = float(mfydip[2])
Emfy = get_energy(mfyfield)
# molecule in -fz field
mfzfield ="mfz"
mfzdip =[]
mfz_field = process_file(mfzfield,mfzdip)
uxmfz = float(mfzdip[0])
uymfz = float(mfzdip[1])
uzmfz = float(mfzdip[2])
Emfz = get_energy(mfzfield)
# molecule in 2fx field
tfield = "tfx"
tfxdip = []
tfx_field = process_file(tfield,tfxdip)
uxtfx = float(tfxdip[0])
uytfx = float(tfxdip[1])
uztfx = float(tfxdip[2])
Etfx = get_energy(tfield)
# molecule in 2fy field
tfyfield ="tfy"
tfydip =[]
tfy_field = process_file(tfyfield,tfydip)
uxtfy = float(tfydip[0])
uytfy = float(tfydip[1])
uztfy = float(tfydip[2])
Etfy = get_energy(tfyfield)
# molecule in 2fz field
tfzfield ="tfz"
tfzdip =[]
tfz_field = process_file(tfzfield,tfzdip)
uxtfz = float(tfzdip[0])
uytfz = float(tfzdip[1])
uztfz = float(tfzdip[2])
Etfz = get_energy(tfzfield)
# molecule in -2fx field
mtfxfield ="mtfx"
mtfxdip =[]
mtfx_field = process_file(mtfxfield,mtfxdip)
uxmtfx = float(mtfxdip[0])
uymtfx = float(mtfxdip[1])
uzmtfx = float(mtfxdip[2])
Emtfx = get_energy(mtfxfield)
# molecule in -2fy field
mtfyfield ="mtfy"
mtfydip =[]
mtfy_field = process_file(mtfyfield,mtfydip)
uxmtfy = float(mtfydip[0])
uymtfy = float(mtfydip[1])
uzmtfy = float(mtfydip[2])
Emtfy = get_energy(mtfyfield)
# molecule in -2fz field
mtfzfield ="mtfz"
mtfzdip =[]
mtfz_field = process_file(mtfzfield,mtfzdip)
uxmtfz = float(mtfzdip[0])
uymtfz = float(mtfzdip[1])
uzmtfz = float(mtfzdip[2])
Emtfz = get_energy(mtfzfield)
#MOLECULE IN fx,fy field
fxfyf = "fxfy"
fxfydip = []
fxfy_field = process_file(fxfyf,fxfydip)
uxfxfy = float(fxfydip[0])
uyfxfy = float(fxfydip[1])
uzfxfy = float(fxfydip[2])
Efxfy = get_energy(fxfyf)
uxfyfx, uyfyfx, uzfyfx, Efyfx =  uxfxfy, uyfxfy, uzfxfy, Efxfy
#molecule in -fx,fy field
mfxfyf = "mfxfy"
mfxfydip = []
mfxfy_field = process_file(mfxfyf,mfxfydip)
uxmfxfy = float(mfxfydip[0])
uymfxfy = float(mfxfydip[1])
uzmfxfy = float(mfxfydip[2])
Emfxfy = get_energy(mfxfyf)
uxfymfx, uyfymfx, uzfymfx, Efymfx = uxmfxfy, uymfxfy, uzmfxfy, Emfxfy
#molecule in fx,-fy field
fxmfyf = "fxmfy"
fxmfydip = []
fxmfy_field = process_file(fxmfyf,fxmfydip)
uxfxmfy = float(fxmfydip[0])
uyfxmfy = float(fxmfydip[1])
uzfxmfy = float(fxmfydip[2])
Efxmfy = get_energy(fxmfyf)
uxmfyfx, uymfyfx, uzmfyfx, Emfyfx = uxfxmfy, uyfxmfy, uzfxmfy, Efxmfy
#molecule in -fx,-fy field
mfxmfyf = "mfxmfy"
mfxmfydip = []
mfxmfy_field = process_file(mfxmfyf,mfxmfydip)
uxmfxmfy = float(mfxmfydip[0])
uymfxmfy = float(mfxmfydip[1])
uzmfxmfy = float(mfxmfydip[2])
Emfxmfy = get_energy(mfxmfyf)
uxmfymfx, uymfymfx, uzmfymfx, Emfymfx = uxmfxmfy, uymfxmfy, uzmfxmfy, Emfxmfy
#MOLECULE IN (fx,fz) field
fxfzf = "fxfz"
fxfzdip = []
fxfz_field = process_file(fxfzf,fxfzdip)
uxfxfz = float(fxfzdip[0])
uyfxfz = float(fxfzdip[1])
uzfxfz = float(fxfzdip[2])
Efxfz = get_energy(fxfzf)
uxfzfx, uyfzfx, uzfzfx, Efzfx = uxfxfz, uyfxfz, uzfxfz, Efxfz
#molecule in (-fx,fz) field
mfxfzf = "mfxfz"
mfxfzdip = []
mfxfz_field = process_file(mfxfzf,mfxfzdip)
uxmfxfz = float(mfxfzdip[0])
uymfxfz = float(mfxfzdip[1])
uzmfxfz = float(mfxfzdip[2])
Emfxfz = get_energy(mfxfzf)
uxfzmfx, uyfzmfx, uzfzmfx, Efzmfx = uxmfxfz, uymfxfz, uzmfxfz, Emfxfz
#molecule in fx,-fz field
fxmfzf = "fxmfz"
fxmfzdip = []
fxmfz_field = process_file(fxmfzf,fxmfzdip)
uxfxmfz = float(fxmfzdip[0])
uyfxmfz = float(fxmfzdip[1])
uzfxmfz = float(fxmfzdip[2])
Efxmfz = get_energy(fxmfzf)
uxmfzfx, uymfzfx, uzmfzfx, Emfzfx = uxfxmfz, uyfxmfz, uzfxmfz, Efxmfz
#molecule in -fx,-fz field
mfxmfzf = "mfxmfz"
mfxmfzdip = []
mfxmfz_field = process_file(mfxmfzf,mfxmfzdip)
uxmfxmfz = float(mfxmfzdip[0])
uymfxmfz = float(mfxmfzdip[1])
uzmfxmfz = float(mfxmfzdip[2])
Emfxmfz = get_energy(mfxmfzf)
uxmfzmfx, uymfzmfx, uzmfzmfx, Emfzmfx = uxmfxmfz, uymfxmfz, uzmfxmfz, Emfxmfz
#molecule fy,fz field
fyfzf = "fyfz"
fyfzdip = []
fyfz_field = process_file(fyfzf,fyfzdip)
uxfyfz = float(fyfzdip[0])
uyfyfz = float(fyfzdip[1])
uzfyfz = float(fyfzdip[2])
Efyfz = get_energy(fyfzf)
uxfzfy, uyfzfy, uzfzfy, Efzfy = uxfyfz, uyfyfz, uzfyfz, Efyfz
#molecule -fy,fz field
mfyfzf = "mfyfz"
mfyfzdip = []
mfyfz_field = process_file(mfyfzf,mfyfzdip)
uxmfyfz = float(mfyfzdip[0])
uymfyfz = float(mfyfzdip[1])
uzmfyfz = float(mfyfzdip[2])
Emfyfz = get_energy(mfyfzf)
uxfzmfy, uyfzmfy, uzfzmfy, Efzmfy = uxmfyfz, uymfyfz, uzmfyfz, Emfyfz
#molecule fy,-fz field
fymfzf = "fymfz"
fymfzdip = []
fymfz_field = process_file(fymfzf,fymfzdip)
uxfymfz = float(fymfzdip[0])
uyfymfz = float(fymfzdip[1])
uzfymfz = float(fymfzdip[2])
Efymfz = get_energy(fymfzf)
uxmfzfy, uymfzfy, uzmfzfy, Emfzfy = uxfymfz, uyfymfz, uzfymfz, Efymfz
#molecule -fy,-fz field
mfymfzf = "mfymfz"
mfymfzdip = []
mfymfz_field = process_file(mfymfzf,mfymfzdip)
uxmfymfz = float(mfymfzdip[0])
uymfymfz = float(mfymfzdip[1])
uzmfymfz = float(mfymfzdip[2])
Emfymfz = get_energy(mfymfzf)
uxmfzmfy, uymfzmfy, uzmfzmfy, Emfzmfy = uxmfymfz, uymfymfz, uzmfymfz, Emfymfz
#molecule (2fy,2fx) field
tfytfxf = "tfytfx"
tfytfxdip = []
tfytfx_field = process_file(tfytfxf,tfytfxdip)
uxtfytfx = float(tfytfxdip[0])
uytfytfx = float(tfytfxdip[1])
uztfytfx = float(tfytfxdip[2])
Etfytfx = get_energy(tfytfxf)
uxtfxtfy, uytfxtfy, uztfxtfy, Etfxtfy = uxtfytfx, uytfytfx, uztfytfx, Etfytfx
#molecule (-2fy,2fx) field
mtfytfxf = "mtfytfx"
mtfytfxdip = []
mtfytfx_field = process_file(mtfytfxf,mtfytfxdip)
uxmtfytfx = float(mtfytfxdip[0])
uymtfytfx = float(mtfytfxdip[1])
uzmtfytfx = float(mtfytfxdip[2])
Emtfytfx = get_energy(mtfytfxf)
uxtfxmtfy, uytfxmtfy, uztfxmtfy, Etfxmtfy = uxmtfytfx, uymtfytfx, uzmtfytfx, Emtfytfx
#molecule (2fy,-2fx) field
tfymtfxf = "tfymtfx"
tfymtfxdip = []
tfymtfx_field = process_file(tfymtfxf,tfymtfxdip)
uxtfymtfx = float(tfymtfxdip[0])
uytfymtfx = float(tfymtfxdip[1])
uztfymtfx = float(tfymtfxdip[2])
Etfymtfx = get_energy(tfymtfxf)
uxmtfxtfy, uymtfxtfy, uzmtfxtfy, Emtfxtfy = uxtfymtfx, uytfymtfx, uztfymtfx, Etfymtfx
#molecule (-2fy,-2fx) field
mtfymtfxf = "mtfymtfx"
mtfymtfxdip = []
mtfymtfx_field = process_file(mtfymtfxf,mtfymtfxdip)
uxmtfymtfx = float(mtfymtfxdip[0])
uymtfymtfx = float(mtfymtfxdip[1])
uzmtfymtfx = float(mtfymtfxdip[2])
Emtfymtfx = get_energy(mtfymtfxf)
uxmtfxmtfy, uymtfxmtfy, uzmtfxmtfy, Emtfxmtfy = uxmtfymtfx, uymtfymtfx, uzmtfymtfx, Emtfymtfx
#molecule (2fz,2fy) field
tfztfyf = "tfztfy"
tfztfydip = []
tfztfy_field = process_file(tfztfyf,tfztfydip)
uxtfztfy = float(tfztfydip[0])
uytfztfy = float(tfztfydip[1])
uztfztfy = float(tfztfydip[2])
Etfztfy = get_energy(tfztfyf)
uxtfytfz, uytfytfz, uztfytfz, Etfytfz = uxtfztfy, uytfztfy, uztfztfy, Etfztfy
#molecule -2fz,2fy field
mtfztfyf = "mtfztfy"
mtfztfydip = []
mtfztfy_field = process_file(mtfztfyf,mtfztfydip)
uxmtfztfy = float(mtfztfydip[0])
uymtfztfy = float(mtfztfydip[1])
uzmtfztfy = float(mtfztfydip[2])
Emtfztfy = get_energy(mtfztfyf)
uxtfymtfz, uytfymtfz, uztfymtfz, Etfymtfz = uxmtfztfy, uymtfztfy, uzmtfztfy, Emtfztfy
#molecule 2fz,-2fy field
tfzmtfyf = "tfzmtfy"
tfzmtfydip = []
tfzmtfy_field = process_file(tfzmtfyf,tfzmtfydip)
uxtfzmtfy = float(tfzmtfydip[0])
uytfzmtfy = float(tfzmtfydip[1])
uztfzmtfy = float(tfzmtfydip[2])
Etfzmtfy = get_energy(tfzmtfyf)
uxmtfytfz, uymtfytfz, uzmtfytfz, Emtfytfz = uxtfzmtfy, uytfzmtfy, uztfzmtfy, Etfzmtfy
#molecule -2fz,-2fy field
mtfzmtfyf = "mtfzmtfy"
mtfzmtfydip = []
mtfzmtfy_field = process_file(mtfzmtfyf,mtfzmtfydip)
uxmtfzmtfy = float(mtfzmtfydip[0])
uymtfzmtfy = float(mtfzmtfydip[1])
uzmtfzmtfy = float(mtfzmtfydip[2])
Emtfzmtfy = get_energy(mtfzmtfyf)
uxmtfymtfz, uymtfymtfz, uzmtfymtfz, Emtfymtfz = uxmtfzmtfy, uymtfzmtfy, uzmtfzmtfy, Emtfzmtfy
#molecule (2fz,2fx) field
tfztfxf = "tfztfx"
tfztfxdip = []
tfztfx_field = process_file(tfztfxf,tfztfxdip)
uxtfztfx = float(tfztfxdip[0])
uytfztfx = float(tfztfxdip[1])
uztfztfx = float(tfztfxdip[2])
Etfztfx = get_energy(tfztfxf)
uxtfxtfz, uytfxtfz, uztfxtfz, Etfxtfz = uxtfztfx, uytfztfx, uztfztfx, Etfztfx
#molecule (-2fz,2fx) field
mtfztfxf = "mtfztfx"
mtfztfxdip = []
mtfztfx_field = process_file(mtfztfxf,mtfztfxdip)
uxmtfztfx = float(mtfztfxdip[0])
uymtfztfx = float(mtfztfxdip[1])
uzmtfztfx = float(mtfztfxdip[2])
Emtfztfx = get_energy(mtfztfxf)
uxtfxmtfz, uytfxmtfz, uztfxmtfz, Etfxmtfz = uxmtfztfx, uymtfztfx, uzmtfztfx, Emtfztfx
#molecule (2fz,-2fx) field
tfzmtfxf = "tfzmtfx"
tfzmtfxdip = []
tfzmtfx_field = process_file(tfzmtfxf,tfzmtfxdip)
uxtfzmtfx = float(tfzmtfxdip[0])
uytfzmtfx = float(tfzmtfxdip[1])
uztfzmtfx = float(tfzmtfxdip[2])
Etfzmtfx = get_energy(tfzmtfxf)
uxmtfxtfz, uymtfxtfz, uzmtfxtfz, Emtfxtfz = uxtfzmtfx, uytfzmtfx, uztfzmtfx, Etfzmtfx
#molecule (-2fz,-2fx) field
mtfzmtfxf = "mtfzmtfx"
mtfzmtfxdip = []
mtfzmtfx_field = process_file(mtfzmtfxf,mtfzmtfxdip)
uxmtfzmtfx = float(mtfzmtfxdip[0])
uymtfzmtfx = float(mtfzmtfxdip[1])
uzmtfzmtfx = float(mtfzmtfxdip[2])
Emtfzmtfx = get_energy(mtfzmtfxf)
uxmtfxmtfz, uymtfxmtfz, uzmtfxmtfz, Emtfxmtfz = uxmtfzmtfx, uymtfzmtfx, uzmtfzmtfx, Emtfzmtfx
#no field
nfield = " "
nfdip = []
nf_field = process_file(nfield,nfdip)
uxnf = float(nfdip[0])
uynf = float(nfdip[1])
uznf = float(nfdip[2])
E = get_energy(nfield)
#get the electric field
F = get_field(f)
Fx = F*0.707106781#0.5 #F*sin(45)*cos(45)
Fy = F*0.707106781#0.5 #F*sin(45)*cos(45)
Fz = F*0.707106781   #F*cos(45)
rFxz, rFyz = 1.0/Fz,1.0/Fz
rFxFys = 1.0/(Fx*Fy*Fy)
rFyFxs = rFxFys
rFxFzs = 1.0/(Fz*Fz*Fz)
rFyFzs = rFxFzs
rFzFxs = 1.0/(Fz*Fz*Fz)
rFzFys = rFzFxs
rFxFx =1.0/(Fx*Fx)
rFyFy =rFxFx
rFzFz = 1.0/(Fz*Fz)
rfx, rfy, rfz = 1.0/Fx, 1.0/Fy, 1.0/Fz
output = open("staticdir/static_result.mad","w")
output.write("\n\n              Dipole Moment Method                   \n\n")
output.write("Field intensity = %16.8f   \n"%F)
#compute the components of polarizability
tt, ot =  2.0/3.0, 1.0/12.0
axx = rfx*(tt*(uxfx - uxmfx) - ot*(uxtfx - uxmtfx))
ayy = rfy*(tt*(uyfy - uymfy) - ot*(uytfy - uymtfy))
azz = rfz*(tt*(uzfz - uzmfz) - ot*(uztfz - uzmtfz))
am = (1.0/3.0)*(axx + ayy + azz)
output.write("\n\n Diagonal Components of Polarizability \n\n     am = mean polarizability\n")
output.write("axx = %16.8f\nayy = %16.8f\nazz = %16.8f\nam = %16.8f"%(axx,ayy,azz,am))
#compute the components of dipole
rs = 1.0/6.0
ux = (tt*(uxfx + uxmfx) - rs*(uxtfx + uxmtfx))
uy = (tt*(uyfy + uymfy) - rs*(uytfy + uymtfy))
uz = (tt*(uzfz + uzmfz) - rs*(uztfz + uzmtfz))
output.write("\n\n       Dipole Moment          \n\n  ")
output.write("ux = %16.8f \nuy = %16.8f \nuz = %16.8f"%(ux,uy,uz))
#compute off-diagonal components of polarizability
axy = rfy*(tt*(uxfy - uxmfy) - ot*(uxtfy - uxmtfy))
ayx = rfx*(tt*(uyfx - uymfx) - ot*(uytfx - uymtfx))
azx = rfz*(tt*(uzfx - uzmfx) - ot*(uztfx - uzmtfx))
axz = rfz*(tt*(uxfz - uxmfz) - ot*(uxtfz - uxmtfz))
ayz = rfz*(tt*(uyfz - uymfz) - ot*(uytfz - uymtfz))
azy = rfz*(tt*(uzfy - uzmfy) - ot*(uztfy - uzmtfy))
output.write("\n\n     off-Diagonal Components of Polarizability          \n\n")
output.write("axy = %16.8f\nayx = %16.8f\naxz = %16.8f\nayz = %16.8f\nazx = %16.8f\nazy = %16.8f "%(axy,ayx,axz,ayz,azx,azy))
#compute the diagonal components of the hyper-polarizability
th = 1.0/3.0
Bxxx =rfx*rfx*th*(uxtfx + uxmtfx - uxfx -uxmfx)
Byyy =rfy*rfy*th*(uytfy + uymtfy - uyfy -uymfy)
Bzzz =rfz*rfz*th*(uztfz + uzmtfz -uzfz -uzmfz)
output.write("\n\n     Diagonal Components of hyperPolarizability     \n\n")
output.write("Bxxx = %16.8f\nByyy = %16.8f\nBzzz = %16.8f"%(Bxxx,Byyy,Bzzz))
#compute the off-diagonal components of the hyper-polarizability
Bxzz =rFzFz*th*(uxtfz + uxmtfz - uxfz -uxmfz)
Bxyy =rFyFy*th*(uxtfy + uxmtfy - uxfy -uxmfy)
Byzz =rFzFz*th*(uytfz + uymtfz -uyfz -uymfz)
Byxx =rFxFx*th*(uytfx + uymtfx - uyfx -uymfx)
Bzyy =rFzFz*th*(uztfy + uzmtfy - uzfy -uzmfy)
Bzxx =rFzFz*th*(uztfx + uzmtfx -uzfx -uzmfx)
output.write("\n\n     off-Diagonal Components of hyperPolarizability          \n\n")
output.write("Bxyy = %16.8f\nByxx = %16.8f\nBxzz = %16.8f\nByzz = %16.8f\nBzxx = %16.8f\nBzyy = %16.8f"%(Bxyy,Byxx,Bxzz,Byzz,Bzxx,Bzyy))
#compute the scalar projection of hyperpolarizability onto the dipole moment basis
Bx= Bxxx + Bxyy + Bxzz
By= Byxx + Byyy + Byzz
Bz= Bzxx + Bzyy + Bzzz
Bu = Bx*ux + By*uy + Bz*uz
unorm = math.sqrt(ux*ux + uy*uy + uz*uz)
cof =(3.0/5.0)
if unorm == 0.0:
    Bmu = 0.0
else:
    Bmu = cof*(Bu/unorm)
output.write("\n\n     catesian components of hyperpolarizability    \n\n")
output.write("Bx = %16.8f\nBy = %16.8f\nBz = %16.8f\n\ndipole norm = %16.8f\nscalar projection of B = %16.8f "%(Bx,By,Bz,unorm,Bmu))
#compute the diagonal components of the hyper-hyperpolarizability
h = 1.0/2.0
Gxxxx =rfx*rfx*rfx*(h*(uxtfx - uxmtfx) - uxfx +uxmfx)
Gyyyy =rfy*rfy*rfy*(h*(uytfy - uymtfy) - uyfy +uymfy)
Gzzzz =rfz*rfz*rfz*(h*(uztfz - uzmtfz) - uzfz +uzmfz)
output.write("\n\n     Diagonal Components of hyper-hyperPolarizability     \n\n")
output.write("Gxxxx = %16.8f\nGyyyy = %16.8f\nGzzzz = %16.8f"%(Gxxxx,Gyyyy,Gzzzz))
#Compute the off-diagonal components of hyper-hyperpolarizability
Gxxyy = rFxFys*(h*(uxfxfy - uxmfxfy + uxfxmfy - uxmfxmfy) - (uxfx - uxmfx))
Gxxzz = rFxFzs*(h*(uxfxfz - uxmfxfz + uxfxmfz - uxmfxmfz) - (uxfx - uxmfx))
Gyyxx = rFyFxs*(h*(uyfyfx - uymfyfx + uyfymfx - uymfymfx) - (uyfy - uymfy))
Gyyzz = rFyFzs*(h*(uyfyfz - uymfyfz + uyfymfz - uymfymfz) - (uyfy - uymfy))
Gzzxx = rFzFxs*(h*(uzfzfx - uzmfzfx + uzfzmfx - uzmfzmfx) - (uzfz - uzmfz))
Gzzyy = rFzFys*(h*(uzfzfy - uzmfzfy + uzfzmfy - uzmfzmfy) - (uzfz - uzmfz))
Gm = (1.0/5.0)*(Gxxxx + Gyyyy + Gzzzz + 2.0*(Gxxyy + Gyyzz + Gzzxx))
output.write("\n\n     off-Diagonal Components of hyper-hyperPolarizability     \n\n")
output.write("Gxxyy = %16.8f\nGxxzz = %16.8f\nGyyxx = %16.8f\nGyyzz = %16.8f\nGzzxx = %16.8f\nGzzyy = %16.8f"%(Gxxyy,Gxxzz,Gyyxx,Gyyzz,Gzzxx,Gzzyy))
output.write("\nGm = %16.8f"%Gm)
#USING ENERGY TO COMPUTE THE OPTICAL PROPERTIES

#compute the components of polarizability
output.write("\n\n               Energy Method                       \n\n")
fte, fthe, ote =  5.0/2.0, 4.0/3.0, 1.0/12.0
axx = rfx*rfx*(fte*E - fthe*(Efx + Emfx) + ote*(Etfx + Emtfx))
ayy = rfy*rfy*(fte*E - fthe*(Efy + Emfy) + ote*(Etfy + Emtfy))
azz = rfz*rfz*(fte*E - fthe*(Efz + Emfz) + ote*(Etfz + Emtfz))
am = (1.0/3.0)*(axx + ayy + azz)
output.write("\n\n Diagonal Components of Polarizability \n\n     am = mean polarizability\n")
output.write("axx = %16.8f\nayy = %16.8f\nazz = %16.8f\nam = %16.8f"%(axx,ayy,azz,am))
#compute the components of dipole
tthe = 2.0/3.0
ux = (tthe*(Efx - Emfx) + ote*(Etfx - Emtfx))
uy = (tthe*(Efy - Emfy) + ote*(Etfy - Emtfy))
uz = (tthe*(Efz - Emfz) + ote*(Etfz - Emtfz))
output.write("\n\n       Dipole Moment          \n\n  ")
output.write("ux = %16.8f \nuy = %16.8f \nuz = %16.8f"%(ux,uy,uz))
#compute off-diagonal components of polarizability
ofe, othe = 1.0/48.0, 1.0/3.0
axy = rfx*rfy*(ofe*(Etfxtfy - Etfxmtfy - Emtfxtfy + Emtfxmtfy) - othe*(Efxfy - Efxmfy - Emfxfy + Emfxmfy))
ayx = rfy*rfx*(ofe*(Etfytfx - Etfymtfx - Emtfytfx + Emtfymtfx) - othe*(Efyfx - Efymfx - Emfyfx + Emfymfx))
azx = rfz*rfz*(ofe*(Etfztfx - Etfzmtfx - Emtfztfx + Emtfzmtfx) - othe*(Efzfx - Efzmfx - Emfzfx + Emfzmfx))
axz = rfz*rfz*(ofe*(Etfxtfz - Etfxmtfz - Emtfxtfz + Emtfxmtfz) - othe*(Efxfz - Efxmfz - Emfxfz + Emfxmfz))
ayz = rfz*rfz*(ofe*(Etfytfz - Etfymtfz - Emtfytfz + Emtfymtfz) - othe*(Efyfz - Efymfz - Emfyfz + Emfymfz))
azy = rfz*rfz*(ofe*(Etfztfy - Etfzmtfy - Emtfztfy + Emtfzmtfy) - othe*(Efzfy - Efzmfy - Emfzfy + Emfzmfy))
output.write("\n\n     off-Diagonal Components of Polarizability          \n\n")
output.write("axy = %16.8f\nayx = %16.8f\naxz = %16.8f\nayz = %16.8f\nazx = %16.8f\nazy = %16.8f "%(axy,ayx,axz,ayz,azx,azy))
#compute the diagonal components of the hyper-polarizability
he = 1.0/2.0
Bxxx =rfx*rfx*rfx*(Efx - Emfx - he*(Etfx -Emtfx))
Byyy =rfy*rfy*rfy*(Efy - Emfy - he*(Etfy -Emtfy))
Bzzz =rfz*rfz*rfz*(Efz - Emfz - he*(Etfz -Emtfz))
output.write("\n\n     Diagonal Components of hyperPolarizability     \n\n")
output.write("Bxxx = %16.8f\nByyy = %16.8f\nBzzz = %16.8f"%(Bxxx,Byyy,Bzzz))
#compute the off-diagonal components of the hyper-polarizability
Bxzz = rfz*rfz*rfz*(he*(Emfxmfz - Efxfz + Emfxfz - Efxmfz) + Efx - Emfx)
Bxyy = rfx*rfy*rfy*(he*(Emfxmfy - Efxfy + Emfxfy - Efxmfy) + Efx - Emfx)
Byzz = rfz*rfz*rfz*(he*(Emfymfz - Efyfz + Emfyfz - Efymfz) + Efy - Emfy)
Byxx = rfy*rfx*rfx*(he*(Emfymfx - Efyfx + Emfyfx - Efymfx) + Efy - Emfy)
Bzyy = rfz*rfz*rfz*(he*(Emfzmfy - Efzfy + Emfzfy - Efzmfy) + Efz - Emfz)
Bzxx = rfz*rfz*rfz*(he*(Emfzmfx - Efzfx + Emfzfx - Efzmfx) + Efz - Emfz)
output.write("\n\n     off-Diagonal Components of hyperPolarizability          \n\n")
output.write("Bxyy = %16.8f\nByxx = %16.8f\nBxzz = %16.8f\nByzz = %16.8f\nBzxx = %16.8f\nBzyy = %16.8f"%(Bxyy,Byxx,Bxzz,Byzz,Bzxx,Bzyy))
#compute the scalar projection of hyperpolarizability onto the dipole moment basis
Bx= Bxxx + Bxyy + Bxzz
By= Byxx + Byyy + Byzz
Bz= Bzxx + Bzyy + Bzzz
Bu = Bx*ux + By*uy + Bz*uz
unorm = math.sqrt(ux*ux + uy*uy + uz*uz)
cof =(3.0/5.0)
if unorm == 0.0:
    Bmu = 0.0
else:
    Bmu = cof*(Bu/unorm)
output.write("\n\n     catesian components of hyperpolarizability    \n\n")
output.write("Bx = %16.8f\nBy = %16.8f\nBz = %16.8f\n\ndipole norm = %16.8f\nscalar projection of B = %16.8f "%(Bx,By,Bz,unorm,Bmu))
#compute the diagonal components of the hyper-hyperpolarizability
Gxxxx = rfx*rfx*rfx*rfx*(4.0*(Efx + Emfx) - 6.0*E - Etfx - Emtfx)
Gyyyy = rfy*rfy*rfy*rfy*(4.0*(Efy + Emfy) - 6.0*E - Etfy - Emtfy)
Gzzzz = rfz*rfz*rfz*rfz*(4.0*(Efz + Emfz) - 6.0*E - Etfz - Emtfz)
output.write("\n\n     Diagonal Components of hyper-hyperPolarizability     \n\n")
output.write("Gxxxx = %16.8f\nGyyyy = %16.8f\nGzzzz = %16.8f"%(Gxxxx,Gyyyy,Gzzzz))
#Compute the off-diagonal components of hyper-hyperpolarizability
Gxxyy = rFxFys*rfx*(2.0*(Efx + Emfx + Efy + Emfy) - 4.0*E - Efxfy - Emfxfy - Efxmfy - Emfxmfy)
Gxxzz = rFxFzs*rfz*(2.0*(Efx + Emfx + Efz + Emfz) - 4.0*E - Efxfz - Emfxfz - Efxmfz - Emfxmfz)
Gyyxx = rFyFxs*rfy*(2.0*(Efy + Emfy + Efx + Emfx) - 4.0*E - Efyfx - Emfyfx - Efymfx - Emfymfx)
Gyyzz = rFyFzs*rfz*(2.0*(Efy + Emfy + Efz + Emfz) - 4.0*E - Efyfz - Emfyfz - Efymfz - Emfymfz)
Gzzxx = rFzFxs*rfz*(2.0*(Efz + Emfz + Efx + Emfx) - 4.0*E - Efzfx - Emfzfy - Efzmfx - Emfzmfx)
Gzzyy = rFzFys*rfz*(2.0*(Efz + Emfz + Efy + Emfy) - 4.0*E - Efzfy - Emfzfy - Efzmfy - Emfzmfy)
Gm = (1.0/5.0)*(Gxxxx + Gyyyy + Gzzzz + 2.0*(Gxxyy + Gyyzz + Gzzxx))
output.write("\n\n     off-Diagonal Components of hyper-hyperPolarizability     \n\n")
output.write("Gxxyy = %16.8f\nGxxzz = %16.8f\nGyyxx = %16.8f\nGyyzz = %16.8f\nGzzxx = %16.8f\nGzzyy = %16.8f"%(Gxxyy,Gxxzz,Gyyxx,Gyyzz,Gzzxx,Gzzyy))
output.write("\nGm = %16.8f"%Gm)
output.write(" \n Debugging \n Efx = %10.6f Efy = %10.6f  Efz = %10.6f\n Emfx = %10.6f Emfy = %10.6f  Emfz = %10.6f\n Emtfx = %10.6f Emtfy = %10.6f  Emtfz = %10.6f\n"%(Efx, Efy,Efz,Emfx, Emfy,Emfz,Emtfx, Emtfy,Emtfz))
output.write(" \n Debugging \n uxfx = %10.6f uyfy = %10.6f  uzfz = %10.6f\n uxmfx = %10.6f uymfy = %10.6f  uzmfz = %10.6f\n uxmtfx = %10.6f uymtfy = %10.6f  uzmtfz = %10.6f\n"%(uxfx, uyfy,uzfz,uxmfx, uymfy,uzmfz,uxmtfx, uymtfy,uzmtfz))
output.write(" \n Debugging \n uyfyfx = %10.6f\n uyfymfx = %10.6f uymfyfx = %10.6f  uymfymfx = %10.6f\n"%(uyfyfx, uyfymfx,uymfyfx,uymfymfx))

output.close()


