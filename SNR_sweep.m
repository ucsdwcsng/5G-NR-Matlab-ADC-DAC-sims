input_SNRs = 5:5:40;
num_seeds = 5;
interp_fac = 0.75;

out_snr_mat = zeros([numel(input_SNRs),num_seeds]);
ref_snr_mat = zeros([numel(input_SNRs),num_seeds]);
out_snr_mat_diff = zeros([numel(input_SNRs),num_seeds]);
ref_snr_mat_diff = zeros([numel(input_SNRs),num_seeds]);

for i=1:1:numel(input_SNRs)
    for j=1:1:num_seeds
        disp("Starting SNR sim: "+num2str(input_SNRs(i))+" dB, and seed: "+j)
        [out_snr_mat(i,j),out_snr_mat_diff(i,j),diag_snr_mat(i,j),diag_snr_mat_diff(i,j),...
            ref_snr_mat(i,j),ref_snr_mat_diff(i,j)]= output_snr_recovered(input_SNRs(i),j,interp_fac);
    end
end
%%
mean_across_seed_outsnr = mean(out_snr_mat,2);
mean_across_seed_refsnr = mean(ref_snr_mat,2);
mean_across_seed_diagsnr = mean(diag_snr_mat,2);


plot(input_SNRs,mean_across_seed_outsnr,'-','LineWidth',3);
plot(input_SNRs,mean_across_seed_diagsnr,'--','LineWidth',3);
hold on
plot(input_SNRs,mean_across_seed_refsnr,'--','LineWidth',3);
plot(input_SNRs,input_SNRs,':','LineWidth',2);
xlabel("Input SNR")
ylabel("Output SNR")
grid on
grid minor
legend(["Recon SNR (ADC4x)","Recon SNR (ADC4x, Diag Eq)","DBF SNR (ADC1x)","Ideal"])
save("./Data/snr_recon_output_interpfac_"+num2str(interp_fac)+"_numseeds_"+num2str(num_seeds)+".mat",'out_snr_mat','ref_snr_mat','diag_snr_mat')