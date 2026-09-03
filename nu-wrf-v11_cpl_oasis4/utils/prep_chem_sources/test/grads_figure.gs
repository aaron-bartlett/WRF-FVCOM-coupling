*'reinit'
*'open Queima_tracer-T-2006-08-10-000000-g1.ctl'
*'open Queima_tracer-ll-T-2005-06-27-000000-g1.ctl'
*'open Queima_tracer-ll-0.1dg-T-2005-06-27-000000-g1.ctl'
*'reset'
'set mpdset mres'
'set map 1 1 5'
'set xlopts 1 6 0.15'
'set ylopts 1 6 0.15'
'set  clopts -1 -1 0.12'
dummy = colors()

* grid *************************
*'set lon -180 180'
'c'

'set gxout grfill'
'set parea 1 10 1 8'

'q time'
result=subwrd(result,3)
horai=substr(result,1,3)
diai=substr(result,4,2)
mesi=substr(result,6,3)
anoi=substr(result,9,4)

*'set z 'zlev
'set z 1'

'set grads off'
'set grid off'
'set ylab on'
'set xlint 30'
'set ylint 30'

'set ccols 0  22 44   74  31 45  61 80 72 67 62 40   35'
'set clevs  0 0.01 0.05 0.1 0.25 0.5 1  10  50 100 500 750 1000 1250'
*'set clevs  0 0.01 0.05 0.1 0.25 0.5 1  10  50 100 500 1000'
*'set clevs  0 0.01 0.025  0.05 0.1 0.25 0.5 1  1.5 2.5 5 10 '
*'d co_bburn.2*1.e5'
'd co_antro.1*1.e6'
'run cbar'

*'draw title RETRO Anthropogenic global emission \ carbon monoxide (1e-7 kg/m^2, 0.5x0.5 degree) @ 10AUG2006'

*'draw string 2 7.7 'horai''diai''mesi''anoi' - lev='level/1000.' km''
return
*'printim  'var'_z'zlev'_g'grade'_t'it'.gif x550 y450'
file1='co-global-emission-bb'
'enable print 'file1'.gmf'
'print'
'disable print'
'!gxps -c -i 'file1'.gmf -o 'file1'.ps'
'!convert -rotate 90 'file1'.ps 'file1'.gif'

'!gxeps -Rc -i 'file1'.gmf -o 'file1'.eps'






*'q pos'

return

*************************************
function colors
*karla

'set rgb  22    0    0   74'
'set rgb  30    0    0  122'
'set rgb  38    0    0  170'
'set rgb  60   78   91  255'
'set rgb  74  162  189  255'
'set rgb  44    0    0  206'


'set rgb  31    0   96    0'
'set rgb  36    0  126    0'
'set rgb  45    0  174    0'
'set rgb  50   18  204   21'
'set rgb  61   78  244   91'

'set rgb  35  134    0    0'
'set rgb  40  180    0    0'
'set rgb  47  206    0    0'
'set rgb  52  226   40   15'
'set rgb  57  246   80   30'
'set rgb  62  255  120   45'
'set rgb  67  255  160   60'
'set rgb  72  255  200   75'
'set rgb  80  255  255   99'

return
