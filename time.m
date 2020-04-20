for i =1:10
    % P = radon(phantom(512));
    P = paradon(phantom(512),linspace(0,360-360/256,256),600,1);
end