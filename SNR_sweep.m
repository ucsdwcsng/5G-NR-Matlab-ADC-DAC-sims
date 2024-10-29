input_SNRs = 5:5:40;
num_seeds = 40;
interp_fac = 0.1165;

out_snr_mat = zeros([numel(input_SNRs),num_seeds]);
ref_snr_mat = zeros([numel(input_SNRs),num_seeds]);

for i=1:1:numel(input_SNRs)
    for j=1:1:num_seeds
        disp("Starting SNR sim: "+num2str(input_SNRs(i))+" dB, and seed: "+j)
        [out_snr_mat(i,j),ref_snr_mat(i,j)] = output_snr_recovered(input_SNRs(i),j,interp_fac);
    end
end
%%
mean_across_seed_outsnr = mean(out_snr_mat,2);
mean_across_seed_refsnr = mean(ref_snr_mat,2);

plot(input_SNRs,mean_across_seed_outsnr,'-','LineWidth',3);
hold on
plot(input_SNRs,mean_across_seed_refsnr,'--','LineWidth',3);
plot(input_SNRs,input_SNRs,':','LineWidth',2);
grid on
grid minor
legend(["Recon SNR (ADC4x)","DBF SNR (ADC1x)","Ideal"])
save("./Data/snr_recon_output_interpfac_"+num2str(interp_fac)+"_numseeds_"+num2str(num_seeds)+".mat",'out_snr_mat','ref_snr_mat')