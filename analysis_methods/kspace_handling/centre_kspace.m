function kspace_centred = centre_kspace(kspace)

dims = size(kspace);


kspace = reshape(kspace,dims(1),dims(2),[]);
kspace_centred = zeros(size(kspace));
temp = abs(kspace(:,:,1));

            [~,ind] = (max(temp(:)));
            [I1, I2] = ind2sub( size(temp),ind);           
for n=1:size(kspace,3);
    kspace_centred(:,:,n) = circshift(kspace(:,:,n),[dims(1)/2-I1,dims(2)/2-I2]);
end 

kspace_centred = reshape(kspace_centred,dims);