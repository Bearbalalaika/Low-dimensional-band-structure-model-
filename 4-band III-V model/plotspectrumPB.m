[I, f]= pade_baker_dft_turbo(IO.current.probe(5000:30000, 1, 3), 1/0.0286, 0.5, 0.44, 0.968, 6000, 5);
figure; plot(484./f, abs(I).^2);