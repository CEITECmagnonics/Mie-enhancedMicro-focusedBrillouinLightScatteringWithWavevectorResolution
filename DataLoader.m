function [Output,f]=DataLoader(FileName)
load(sprintf('%s.mat',FileName),'Data');


M=length(Data);
Output=cell(M,1);
for kk=1:M
    
    f=Data{kk}.ENiFe.f;
    Nf=length(f);
    Output{kk}.f=f;
    Output{kk}.a=Data{kk}.a; % parameter 'a' for stripe

    xi=Data{kk}.ENiFe.x;
    yi=Data{kk}.ENiFe.y;
    zi=Data{kk}.ENiFe.z;

    Output{kk}.xi=xi;
    Output{kk}.yi=yi;
    Output{kk}.zi=zi;
    Output{kk}.Ei=reshape(squeeze(Data{kk}.ENiFe.E),[length(xi),length(yi),3,Nf]);

end














    
    




