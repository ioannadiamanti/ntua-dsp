function b = hz2bark(f)
    b = 13*atan(.00076.*f)+3.5*atan((f./7500).^2);
end