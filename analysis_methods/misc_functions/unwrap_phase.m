function unwrapped_phase = unwrap_phase(compleximages)

dims = size(compleximages);
temp = reshape(compleximages,dims(1),dims(2),[]);
unwrapped_phase = zeros(size(temp));
for n=1:size(temp,3)
    [th4noise] = NoiseEstimation(temp(:,:,n));
    % prepare the mask
    th4unwrap = 4.0; % the pixels included in the unwrapping mask
    th4supp = 4.0; % the pixels included in the support mask
    th4stack = 15.0; % the pixels used for stacking the multiple slices
    [mask4unwrap, mask4supp, mask4stack] = mask_generation(temp(:,:,n),th4noise,th4unwrap,th4supp,th4stack);
    %---------------------------------------------------------------------
    [  unwrapped_phase(:,:,n) ] = PUROR2D( temp(:,:,n),mask4unwrap, mask4supp, mask4stack, 0 );
end
unwrapped_phase = reshape(unwrapped_phase,dims);