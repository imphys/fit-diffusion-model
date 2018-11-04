function save_niiTgz( fn, img, T , tempdir)
%    save_niiTgz( fn, nii, [], tempdir);
%    save_niiTgz( fn, img, T, tempdir);
% Saves img to fn
% Convenience interface to immediately store the nifti. Calls make_niiT( img, T)
% Call gz on the stored nifti and remove the nii file.
% Don't include the .gz extension in fn; its generated automatically.
%
% Created by Dirk Poot, Erasmus MC. 19-7-2012

if ~isstruct(img)
    nii = make_niiT(img,  T);
else
    nii = img;
end;

if nargin>3
    [outdir, outfile, outext ] = fileparts(fn);
    tempfile = fullfile( tempdir, [outfile outext]); 
    fprintf('saving: %s', tempfile );
    save_nii(nii, tempfile);
    fprintf('.gz');
    if 0
        % use gzip store in different directory argument
        gzip( tempfile, outdir );
    else
        % zip in place and move (seems to be MUCH faster over network)
        gzip( tempfile );
        fprintf('\nMove gz file to %s', outdir);
        movefile( [tempfile '.gz'] , outdir);
    end;
    delete( tempfile );
else
    fprintf('saving: %s', fn );
    save_nii(nii, fn);
    fprintf('.gz');
    gzip( fn );
    delete( fn );
end;
fprintf('\n');
