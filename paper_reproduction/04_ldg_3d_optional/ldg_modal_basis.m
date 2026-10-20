function phi=ldg_modal_basis(xyz)
phi=zeros(8,size(xyz,1));
for c=0:1
    for b=0:1
        for a=0:1
            ids=[a b c]; v=ones(size(xyz,1),1);
            for d=1:3
                if ids(d)==0, v=v/sqrt(2);
                else, v=v.*(sqrt(3/2)*(2*xyz(:,d)-1)); end
            end
            phi(1+a+2*b+4*c,:)=v';
        end
    end
end
end
