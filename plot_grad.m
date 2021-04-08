function plot_grad(grad, dimension)
    if(strcmpi('x', dimension))
        time = cumsum(arrayfun(@(x) x.time, grad));
        plot(time, arrayfun(@(x) real(x.G), grad)) 
    else
        time = cumsum(arrayfun(@(x) x.time, grad));
        plot(time, arrayfun(@(x) imag(x.G), grad))
    end
end