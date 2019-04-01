function f=get_resp_phase(trial_file,respiration_file)

    %quick function to get respiration phase regressor from a Hilbert transform
    
    X = csvread(trial_file);

    delimiterIn = '\n';
    headerlinesIn = 0;
    phys_data = importdata(respiration_file,delimiterIn,headerlinesIn);
    %remove NaN values (perhaps due to labels)
    phys_data(~any(~isnan(phys_data), 2),:)=[];
    limit_physio_in_seconds=floor((length(phys_data)/50)-4);
    
    
    indices = find(X(:,1,1)>limit_physio_in_seconds);
    X(indices) = NaN;
    X(sum(isnan(X), 2) == 1, :) = [];

      
    
    h = fspecial('gaussian', 25, 10);
    %J = real(ifft2(fft2(phys_data).*fft2(h, size(phys_data),1)));
    J = imfilter(phys_data, h);

    dX=J(3:end)-J(1:end-2);
    xx=zscore(J(2:end-1))-i*zscore(dX);
    phase_method2=[nan;angle(hilbert(xx));nan]; 
 
 
    %get first column
    X=X';
    A=(X(1,:,:))';

    A=ceil(50.*A);
    resulting_inst_phase=phase_method2(A(:));
    

    [pathstr,name,ext] = fileparts(trial_file);
    ss=strcat(pathstr,'/inst_phase.txt');
    fid=fopen(ss,'w');
    fprintf(fid, '%f\n', resulting_inst_phase);
    fclose(fid);
    f=1;

end