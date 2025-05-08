#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TNtuple.h>
#include <iostream>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <fstream>

void read_ntuple_from_file(const char* filename, 
                           const float xmin, 
                           const float xmax, 
                           const float ymin, 
                           const float ymax,
                           const float lumi) {

   /* Open the NTuple file. */
   TFile* file = TFile::Open(filename);
   if (!file || file->IsZombie()) {
      std::cerr << "Error opening file." << std::endl;
      return;
   }
   
   /* Get the NTuple names. */
   TIter next(file->GetListOfKeys());
   TKey* key;
   while ((key = (TKey*)next())) {
        std::cout << "Name: " << key->GetName()
                  << " Class: " << key->GetClassName() << std::endl;
   }

   /* Read the first TTree in the root file. */
   key = (TKey*)file->GetListOfKeys()->At(1);   
   TObject* obj = key->ReadObj();
   /* Check if the key is actually a TTree. */
   if (obj->InheritsFrom("TTree")) {
      TTree* tree = (TTree*)obj;

      /* Print the names of the branches. Use them later in *
       * in the code to assign addresses,                   */
      /*TIter next(tree->GetListOfBranches());
      TBranch* branch;
      while ((branch = (TBranch*)next())) {
         std::cout << "Branch: " << branch->GetName() << std::endl;
      } */

      /* Declare the weight, momentum, position and pdgID. */
      float factor;
      float p[4];
      float x[3];
      float st3_x[3];
      int pdg;
      float nmuons = 0;
      float nantimuons = 0;
   
      float theta_x;
      float theta_y;

      tree->SetBranchAddress("P", &p[0]);
      tree->SetBranchAddress("px", &p[1]);
      tree->SetBranchAddress("py", &p[2]);
      tree->SetBranchAddress("pz", &p[3]);
      tree->SetBranchAddress("x", &x[0]);
      tree->SetBranchAddress("y", &x[1]);
      tree->SetBranchAddress("z", &x[2]);
      tree->SetBranchAddress("st3_x", &st3_x[0]);
      tree->SetBranchAddress("st3_y", &st3_x[1]);
      tree->SetBranchAddress("st3_z", &st3_x[2]);
      tree->SetBranchAddress("factor", &factor);
      tree->SetBranchAddress("pdg", &pdg);

      /* Create the logarithmically spaced bins. */
      float logmin = -3;
      float logmax = 0;
      int binsPerDecade = 7;
      int decades = logmax - logmin;
      int nbins = binsPerDecade * decades;
      float* bins = new float[nbins + 1];
      float binwidth = (float)(logmax - logmin) / (float)nbins;
      for (int i = 0; i <= nbins; i = i + 1) {
         bins[i] = std::pow(10, logmin + i * binwidth);
      }

      /* Book a histogram to visualise the extracted muon flux. */
      TH1F* hmu = new TH1F("muon", "Muon and Anti-Muon Flux at FASER;E_{#mu};N", nbins, bins);
      TH1F* hmub = new TH1F("antimuon", "Anti-muon Flux at FASER;E_{#bar{#mu)};N", nbins, bins);
      //TH1F* h1 = new TH1F("tracks", "Reproduce Figure 84;Momentum [GeV];number of tracks [/2GeV]",100,0,200);
      TH1F* h1 = new TH1F("tracks", "Reproduce Figure 83;Momentum [MeV];events per bin",nbins, bins);

      float target_length = 50;
      float N_avogadro = 6.02e23;
      float rho = 19.3;
      float target_nucleon_density = rho * N_avogadro / 1.0e36;
      float norm = target_nucleon_density * target_length;

      /* Loop over all entries in the NTuple file and write to * 
       * histogram for the flux if passing cuts.               */
      for (int irow = 0; irow < tree->GetEntries(); irow = irow + 1) {
         tree->GetEntry(irow);

         theta_x = std::abs(p[1]/p[3]);
         theta_y = std::abs(p[2]/p[3]);

         //theta_x = std::abs(x[0]/p[2]);
         //theta_y = std::abs(x[1]/p[2]);

         float r;
         /* Fix this to the position of FASERnu in Run 3. */
         float offsetx = 10;
         float offsety = -33;
         r = std::sqrt(std::pow(x[0]-offsetx,2)+std::pow(x[1]-offsety,2));


         if (x[0] < xmax && x[0] > xmin && x[1] < ymax && x[1] > ymin) { 
         //if (x[0] < xmax && x[0] > xmin && x[1] < ymax && x[1] > ymin && theta_x < 0.01 && theta_y < 0.01 && r < 90) {
         //if (r < 90) {
         //if (r < 45.0 ) {
            /*
            std::cout << "factor: " << factor 
                      << " pdg: " << pdg
                      << " E: " << p[0] 
                      << " px: " << p[1]
                      << " py: " << p[2] 
                      << " pz: " << p[3] 
                      << " x: " << x[0] 
                      << " y: " << x[1] 
                      << " z: " << x[2] 
                      << std::endl;
            */

            /* Fill the histograms with the weight factor. Covert MeV to GeV.*/
            if (pdg == 13) {
               hmu->Fill(p[0]/1000/7000, norm * 2/1.3713*factor*lumi);
            //   if (p[0]/1000 > 50 && p[0]/1000 < 1000) nmuons = nmuons + factor;
            }
            if (pdg == -13) {
               hmub->Fill(p[0]/1000/7000, norm * 2/1.3713*factor*lumi);
               if (p[0]/1000 > 50 && p[0]/1000 < 1000) nantimuons = nantimuons + factor;
            }
            nmuons = nmuons + factor;
         }
      }
      printf("Muon 50 GeV < E < 1 TeV : %14.7e\n", nmuons);
      printf("Anti 50 GeV < E < 1 TeV : %14.7e\n", nantimuons);
      std::cout << tree->GetEntries() << std::endl;
      TCanvas* c1 = new TCanvas("c1", "Canvas", 800, 600);
      /* Make the axes logarithmically scaled. */
      c1->SetLogx();
      c1->SetLogy();
      /* Remove the box in the top right corner. */
      gStyle->SetOptStat(0);
      /* Set the colours. */
      hmu->SetLineColor(kBlue);
      hmub->SetLineColor(kRed);
      /* Increase the linewidth. */
      hmu->SetLineWidth(2);
      hmub->SetLineWidth(2);
      h1->SetLineWidth(2);
      /* Draw the lines onto the canvas. */
      hmu->Draw("HIST");
      hmub->Draw("HIST SAME");
      /* Add a legend to the histogram. */
      TLegend* legend = new TLegend(0.15, 0.15, 0.4, 0.3);
      legend->AddEntry(hmu, "Muon flux", "l");
      legend->AddEntry(hmub, "Anti-Muon flux", "l");
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->Draw();
      /* Add a labels for the luminosity */
      TLatex* latex = new TLatex();
      latex->SetNDC();               
      latex->SetTextSize(0.03);      
      TString label = Form("#scale[1.2]{#it{L}_{int}= %d fb^{-1}}", (int)lumi);
      latex->DrawLatex(0.8, 0.95, label);
      c1->SaveAs("muon_flux.pdf");

      std::ofstream lhagridfile("muon_flux_FASERv_Run3_var2_0000.dat");
      if (lhagridfile.is_open()) {
         /* Header of the LHAPDF file. */
         lhagridfile << "PdfType: replica\n";
         lhagridfile << "Format: lhagrid1\n";
         lhagridfile << "---\n";
         /* Write values for xnu here, the central points of the bins. */
         for (int i = 1; i <= hmu->GetNbinsX(); i = i + 1) {
            lhagridfile << hmu->GetBinCenter(i) << " ";             
         }
         lhagridfile << "1.0" << std::endl;
         lhagridfile << "0.1E+001 0.1E+007" << std::endl;
         lhagridfile << "-13 13" << std::endl;
         /* Write the grid, content of the historam. */
         for (int i = 1; i <= hmu->GetNbinsX(); i = i + 1) {
            float valmu = hmu->GetBinContent(i) / hmu->GetBinWidth(i) * hmu->GetBinCenter(i);
            float valmub = hmub->GetBinContent(i) / hmub->GetBinWidth(i) * hmub->GetBinCenter(i);
            lhagridfile << valmub << " " << valmu << std::endl;
            lhagridfile << valmub << " " << valmu << std::endl;
            //std::cout << hmu->GetBinWidth(i) << std::endl;
            //std::cout << hmu->GetBinLowEdge(i) << std::endl;
            std::cout << hmu->GetBinContent(i) << std::endl;

         }
         lhagridfile << "0.0 0.0" << std::endl;
         lhagridfile << "0.0 0.0" << std::endl;
         lhagridfile << "---" << std::endl;
      } else {
         std::cout << "Object is not a TTree." << std::endl;
      }
   }
      
   file->Close();
}

int main() {
   read_ntuple_from_file("../../Downloads/Muon_truth_filtered.root", -115.0, 135, -183, 117, 150.0);
   //read_ntuple_from_file("../../Downloads/Muon_truth_filtered.root", -50, 50, -50, 50, 150.0);
}
