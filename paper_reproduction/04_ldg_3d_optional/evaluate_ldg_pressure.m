function values=evaluate_ldg_pressure(result,xyz)
N=result.N; nd=2*N;
cellidx=min(N,max(1,floor(xyz*N)+1));
ref=xyz*N-(cellidx-1);phi=ldg_nodal_basis(ref);
values=zeros(size(xyz,1),1);
for c=0:1
    for b=0:1
        for a=0:1
            ids=2*cellidx(:,1)-1+a+nd*(2*cellidx(:,2)-2+b)+nd^2*(2*cellidx(:,3)-2+c);
            values=values+result.pressure(ids).*phi(1+a+2*b+4*c,:)';
        end
    end
end
end
