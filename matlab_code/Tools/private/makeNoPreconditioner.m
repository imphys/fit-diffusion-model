function [Rs, mulfun] = makeNoPreconditioner(hessinfo, upperbandw, DM, DG )
Rs.dummy = [];
mulfun = @MulNoPreconditioner;

function [x] = MulNoPreconditioner(x, R)
