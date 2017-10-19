function w = sqrwv ( x, La, Lb )
jdim = size(x,2);
Lh   = (Lb - La)/2;
Lq   = Lh / 2;
Lleft  = La + Lq;
Lright = La + Lh;
w    = zeros(jdim,1);
w    = (x(:)>Lleft & x(:)<Lright).* ( w + 1 );
