ofid
clear fields times AA3
for n=1:50
fields(n) = str2num(object.param.b_evol{1+n});
times(n) = str2num(object.param.t_evol{1+n});
end
n_fields = length(unique(fields));
n_times = object.nmrexperiments/n_fields;
times = times(1:n_fields*n_times);
times = reshape(times,[n_times n_fields]);
n_times =1;
dim1 = object.nmrsamples;
dim2 = object.nmrviews;

AA = reshape(A,[dim1 dim2 2 n_times.*n_fields]);
AA = squeeze(AA(:,:,2,:,:));
AA1 = squeeze(AA(:,:,:));
ofid

AA = reshape(A,[dim1 dim2 2 n_times.*n_fields]);
AA = squeeze(AA(:,:,2,:,:));
AA2 = squeeze(AA(:,:,:));
AA3(:,:,:,1)=AA1;
AA3(:,:,:,2)=AA2;

 windowsize = round(0.35*dim1);
            if mod(windowsize,2)
                windowsize = windowsize+1;
            end

  window1 = kaiser(windowsize,3);
     windowfull = window1*window1';
                        if length(windowfull)>size(AA3,1)
                            winl = length(windowfull);
                            iml = dim1;
                            windowfull =  windowfull(round((winl - iml))/2+1:end - round((winl - iml)/2),round((winl - iml)/2)+1:end - round((winl - iml)/2));
                        end
                        if length(windowfull)<dim1
                            winl = size(windowfull,1);
                            iml = dim1;
                            windowfull = padarray(windowfull,[fix((iml-winl)/2), fix((iml-winl)/2)]);
                        end
                        for i=1:2
                            for j=1:n_times*n_fields
                            AA3(:,:,j,i) =  AA3(:,:,j,i).*windowfull;
                            end
                        end
                        
                        Atemp = reshape(AA3,[dim1 dim2 2*n_times.*n_fields]); 
                        
AAAA(:,:,:) = (iterative_images_correction_v7(Atemp,0.15,20,0.15,[]));

clear Atemp


AAAAA = reshape(AAAA,[dim1 dim2 n_times n_fields 2]);
clear AAAA

for i = 1:n_fields
registeredimages1(:,:,:,i) = register_images(AAAAA(:,:,:,i,1),n_times);
registeredimages2(:,:,:,i) = register_images(AAAAA(:,:,:,i,2),n_times);
end


clear AAAAA AA1 AA2 AA3 
