//
// MIDAS analyzer example 2: ROOT analyzer
//
// K.Olchanski
//

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "v1190unpack.h"

#define DELETE(p) if (p) { delete(p); (p)=NULL; }

class EmmaModule: public TAModuleInterface
{
public:
   void Init(const std::vector<std::string> &args);
   void Finish();
   TARunInterface* NewRun(TARunInfo* runinfo);

   bool fVerboseV1190;
   int fTotalEventCounter;
};

struct EmmaRun: public TARunInterface
{
   EmmaModule* fModule;
   int fCounter;

   // EMMA things go below here:

   int ct;
   TH1D *anode_pathway_multiplicity;

   // A1 - A2
   // [0] = A1B - A2B
   // [1] = A1B - A2M
   // [2] = A1B - A2T
   // [3] = A1M - A2B
   // [4] = A1M - A2M
   // [5] = A1M - A2T
   // [6] = A1T - A2B
   // [7] = A1T - A2M
   // [8] = A1T - A2T
   // [9] = A1 - A2 (All Straight)
   // [10] = A1 - A2 (All Crossing)
   // [11] = A1B - A2B (Copy for plotting purpose)
   // [12] = A1M - A2M (Copy for plotting purpose)
   // [13] = A1T - A2T (Copy for plotting purpose)
   TH1D *a1_a2_diff[14];
   // The X and Y differences
   // [0] = X1R - X1L
   // [1] = Y1B - Y1T
   // [2] = X2R - X2L
   // [3] = Y2B - Y2T
   TH1D *x_y_diff[4];
   // The X and Y sums
   // [0] = X1R + X1L - 1AM
   // [1] = Y1B + Y1T - 1AM
   // [2] = X2R + X2L - 2AM
   // [3] = Y2B + Y2T - 2AM
   TH1D *x_y_sum[4];

   // T2 histograms for diff vs sum
   // TH2F *x_y_diff_vs_sum[4];
   //TH2F *x_y_diff_vs_sum_large[4];
   TH2F *x_y_diff_vs_sum_silly[4];

   // Debugging histogram
   TH1D *m_counts[12];

   TCanvas* fCanvasTdcRaw;
   TCanvas* fCanvasA1A2Summary;
   TCanvas* fCanvasA1A2Detail;
   TCanvas* fCanvasMCounts;
   TCanvas* fCanvasXYDiffs;
   TCanvas* fCanvasXYSums;
   TCanvas* fCanvasXYDiffVsSum;
   TCanvas* fCanvasXYPositions;

   TH1D*    fHTdcTrig;
   TH1D*    fHTdcRaw[64];

   EmmaRun(TARunInfo* runinfo, EmmaModule* m)
      : TARunInterface(runinfo)
   {
      printf("EmmaRun::ctor!\n");
      fModule = m;

      // initialize canvases

      fCanvasTdcRaw = new TCanvas("TDC_raw_data");
      //fCanvasA1A2Summary = new TCanvas("A1_A2_Summary");
      //fCanvasA1A2Detail = new TCanvas("A1_A2_Detail");
      fCanvasMCounts = new TCanvas("M_counts");
      fCanvasXYDiffs = new TCanvas("XY_Diffs");
      fCanvasXYSums = new TCanvas("XY_Sums");
      fCanvasXYDiffVsSum = new TCanvas("XY_Diffs_vs_Sum");
      fCanvasXYPositions = new TCanvas("XY_Positions");

      // initialize histograms

      fHTdcTrig = new TH1D("TDC_trig", "TDC trigger signal", 100, 20000, 21000);

      for (int i=0; i<64; i++) {
         char title[256];
         sprintf(title, "TDC_%d", i);
         fHTdcRaw[i] = new TH1D(title, title, 300, -100, 200);
      }

      // **************************************
      // INITIALIZATION 
      // **************************************
      ct = 0;
      // initialize histograms
      char name[100];
      Float_t xdiffmin=-80;
      Float_t ydiffmin=-150;
      Float_t xdiffrange=450;
      Float_t ydiffrange=310;
      Float_t xsummin=310;
      Float_t ysummin=90;
      Float_t xsumrange=80;
      Float_t ysumrange=xsumrange;
      Float_t diffbin=512;
      Float_t sumbin=100;
      
      // initialize anode_pathway_multiplicity histogram
      anode_pathway_multiplicity = new TH1D("APM","Anode Pathway Multiplicity Index",511,0.5,511.5);
      anode_pathway_multiplicity->SetXTitle("Index");
      
      /*{
         // initialize a1_a2_diff histogram
         const char *title[] = {"A1B - A2B", "A1B - A2M", "A1B - A2T",
                                "A1M - A2B", "A1M - A2M", "A1M - A2T",
                                "A1T - A2B", "A1T - A2M", "A1T - A2T",
                                "A1 - A2 (Straight)", "A1 - A2 (Crossing)",
                                "A1B - A2B", "A1M - A2M", "A1T - A2T"};
         for(int i = 0; i < 14; i++){   
            sprintf(name,"a1_a2_diff_%i",i);
            a1_a2_diff[i] = new TH1D(name,title[i],150,-14.5,15.5);
            a1_a2_diff[i]->SetXTitle("A1-A2 time diff (ns)");
         }
         a1_a2_diff[0]->SetLineColor(2);
         a1_a2_diff[4]->SetLineColor(1);
         a1_a2_diff[8]->SetLineColor(4);
         a1_a2_diff[11]->SetLineColor(2);
         a1_a2_diff[12]->SetLineColor(1);
         a1_a2_diff[13]->SetLineColor(4);
      }*/
      
      {
         // initialize m_count histogram
         const char* title[] = {"Anode Multiplicity", "X1L", "X1R", "Y1B", "Y1T", "Anode & Cathode Measurement Total", "64-Ch Measurement Total"};
         for(int i = 0; i < 7; i++){
            sprintf(name,"m_count_%i",i);
            if(i<5)
               m_counts[i] = new TH1D(name,title[i],5,0,4);
            else
               m_counts[i] = new TH1D(name,title[i],100,0,100);
         }
      }
      
      {
         // initialize x_y_diff histogram
         const char* title[] = {"X1R - X1L", "Y1B - Y1T"};
         for(int i = 0; i < 2; i++){
            sprintf(name,"x_y_diff_%i",i);
            if(i == 0)
               x_y_diff[i] = new TH1D(name,title[i],diffbin,ydiffmin,ydiffmin+ydiffrange);
            if(i == 1)
               x_y_diff[i] = new TH1D(name,title[i],diffbin,xdiffmin,xdiffmin+xdiffrange);
            x_y_diff[i]->SetXTitle("Time Diff (ns)");
         }
      }
      
      {
         // initialize x_y_sum histogram
         const char* title[] = {"X1R + X1L", "Y1T + Y1B"};
         for(int i = 0; i < 2; i++){
            sprintf(name,"x_y_sum_%i",i);
            if(i == 0)
               x_y_sum[i] = new TH1D(name,title[i],sumbin,ysummin,ysummin+ysumrange);
            if(i == 1)
               x_y_sum[i] = new TH1D(name,title[i],sumbin,xsummin,xsummin+xsumrange);
            x_y_sum[i]->SetXTitle("Time Sum (ns)");
         }
      }
      {
         // initialize x_y_diff_vs_sum histogram
         const char*title[] = {"X1R/X1L Diff vs Sum", "Y1B/Y1T Diff vs Sum"};
         for(int i = 0; i < 2; i++){
            sprintf(name,"x_y_diff_vs_sum_silly_%i",i);
            if(i == 0)
               x_y_diff_vs_sum_silly[i] = new TH2F(name,title[i],sumbin,ysummin,ysummin+ysumrange,diffbin,ydiffmin,ydiffmin+ydiffrange);
            if(i == 1)
               x_y_diff_vs_sum_silly[i] = new TH2F(name,title[i],sumbin,xsummin,xsummin+xsumrange,diffbin,xdiffmin,xdiffmin+xdiffrange);
            x_y_diff_vs_sum_silly[i]->SetXTitle("Cathode Time Sum (ns)");      
            x_y_diff_vs_sum_silly[i]->SetYTitle("Cathode Time Differecnce (ns)");
         }
      }
   }

   ~EmmaRun()
   {
      printf("EmmaRun::dtor!\n");

      DELETE(fCanvasTdcRaw);
      DELETE(fCanvasMCounts);
      DELETE(fCanvasXYDiffs);
      DELETE(fCanvasXYSums);
      DELETE(fCanvasXYDiffVsSum);
      DELETE(fCanvasXYPositions);
   }

   void ResetHistograms()
   {
      // ***********************************
      // RESET HISTOGRAMS
      // ***********************************
      for(int i =0; i < 7; i++)
         m_counts[i]->Reset();
      for(int i =0; i < 2; i++) {
         x_y_diff[i]->Reset();
         x_y_sum[i]->Reset();
         x_y_diff_vs_sum_silly[i]->Reset();
      }
   }

   void PlotHistograms(TARunInfo* runinfo)
   {
      printf("PlotHistograms!\n");

      // ********************************
      // DRAW HISTOGRAMS
      // ********************************
      // A1-A2 Summary tab
      //   4 pads:    [3 straight paths overplot]    [straight combined]
      //              [crossing combined]            [anode pathway multiplicity]
      //
      // A1 - A2
      // [0] = A1B - A2B
      // [1] = A1B - A2M
      // [2] = A1B - A2T
      // [3] = A1M - A2B
      // [4] = A1M - A2M
      // [5] = A1M - A2T
      // [6] = A1T - A2B
      // [7] = A1T - A2M
      // [8] = A1T - A2T
      // [9] = A1 - A2 (All Straight)
      // [10] = A1 - A2 (All Crossing)
      // [11] = A1B - A2B (Copy for plotting purpose)
      // [12] = A1M - A2M (Copy for plotting purpose)
      // [13] = A1T - A2T (Copy for plotting purpose)

      if (fCanvasTdcRaw) {
         fCanvasTdcRaw->Clear();
         fCanvasTdcRaw->Divide(2,3);

         int p=0;
         fCanvasTdcRaw->cd(++p);
         fHTdcTrig->Draw();

         for (int i=0; i<4; i++) {
            fCanvasTdcRaw->cd(++p);
            fHTdcRaw[i]->Draw();
         }

         fCanvasTdcRaw->Modified();
         fCanvasTdcRaw->Update();
         printf("update!\n");
      }

      {
         TCanvas* c1 = fCanvasMCounts;
         c1->Clear();
         c1->Divide(4,2);
         c1->cd(1);  m_counts[0]->Draw();
         c1->cd(2);  m_counts[1]->Draw();
         c1->cd(3);  m_counts[6]->Draw();
         c1->cd(4);  m_counts[7]->Draw();
         c1->cd(5);  m_counts[2]->Draw();
         
         c1->Modified();
         c1->Update();
      }
      
      {       
         TCanvas* c1 = fCanvasXYDiffs;
         c1->Clear();
         c1->Divide(2,2);
         for(int i = 0; i < 4; i++){
            c1->cd(1+i);
            x_y_diff[i]->Draw();
         }
         c1->Modified();
         c1->Update();
      }
      
      {       
         TCanvas* c1 = fCanvasXYSums;
         c1->Clear();
         c1->Divide(2,2);
         for(int i = 0; i < 4; i++){
            c1->cd(1+i);
            x_y_sum[i]->Draw();
         }
         c1->Modified();
         c1->Update();
      }
      
      { 
         TCanvas* c1 = fCanvasXYDiffVsSum;
         c1->Clear();
         if(!(c1->GetShowEventStatus()))c1->ToggleEventStatus();
         if(!(c1->GetShowToolBar()))c1->ToggleToolBar();
         c1->Divide(2,2);
         for(int i = 0; i < 4; i++){
            c1->cd(1+i);
            x_y_diff_vs_sum_silly[i]->Draw("colz");
         }
         c1->Modified();
         c1->Update();
      }

      {       
         TCanvas* c1 = fCanvasXYPositions;
         c1->Clear();
         c1->Divide(2,2);
         for(int i = 0; i < 4; i++){
            c1->cd(1+i);
            x_y_diff_vs_sum_silly[i]->Draw("colz");
         }
         c1->Modified();
         c1->Update();
      }

      
   }

   void UpdateHistograms(TARunInfo* runinfo, const v1190event* tdc_data)
   {
      std::vector<int> anode_pathway(9,0);
      std::vector<double> earliest_times(64,999999);
      std::vector<int> counts(64,0);

      int tdc_trig = 0;
      for (unsigned int i=0; i<tdc_data->hits.size(); i++) {
         if (tdc_data->hits[i].trailing) // skip trailing edge hits
            continue;
         int chan = tdc_data->hits[i].channel;
         if (chan != 7) // skip if not TDC trigger channel
            continue;
         tdc_trig = tdc_data->hits[i].measurement;
         break;
      }

      printf("tdc_trig %d\n", tdc_trig);

      fHTdcTrig->Fill(tdc_trig);

      int chan = -1;
      // Seems to be some noise in the measurements.  In the case of multiple
      // measurements for the same channel, get the earliest measurement.
      //    double a1m_earliest = 9999999.0, a2m_earliest = 9999999.0;
      // Vector of the earliest TDC time for each channel
      //    std::vector<double> a1_pulse_time(10,-1.0);
      //    int a1_counter = 0;
      for(unsigned int i = 0; i < tdc_data->hits.size(); i++){ // loop over measurements
         if (tdc_data->hits[i].trailing) // skip trailing edge hits
            continue;
         chan = tdc_data->hits[i].channel;
         double t = (tdc_data->hits[i].measurement-tdc_trig) * 0.100; // convert to nsec
         if (fHTdcRaw[chan]) {
            printf("chan %d, time %f\n", chan, t);
            fHTdcRaw[chan]->Fill(t);
         }
         counts[chan] = counts[chan] + 1;
         if(t < earliest_times[chan])
            earliest_times[chan] = t;
      }

      // ***************************************
      // RAW DATA
      // ***************************************
      // Compute and fill anode pathway multiplicity data
      for (int i=0; i < 3; i++) {
         if(earliest_times[i] < 999990) {
            if(earliest_times[3] < 999990) {
               anode_pathway[0+3*i] = 1;
            }
            if(earliest_times[4] < 999990) {
               anode_pathway[1+3*i] = 1;
            }
            if(earliest_times[5] < 999990) {
               anode_pathway[2+3*i] = 1;
            }
         }
      }
      int apm = 0;
      for (int i=0; i<9; i++)
         apm = apm + anode_pathway[i]*int(pow(2,i));
      // multiplicity of 1 means APM = 1, 2, 4, 8, 16, 32, 64, 128,or 256
      anode_pathway_multiplicity->Fill(apm);
      
      // Fill m_count debugger data
      int multiplicity = 0;
      for (int i=0; i<9; i++)
         multiplicity = multiplicity + anode_pathway[i];
      
      m_counts[0]->Fill(multiplicity);
      m_counts[1]->Fill(counts[4]);
      m_counts[2]->Fill(counts[6]);
      m_counts[3]->Fill(counts[7]);
      m_counts[4]->Fill(counts[8]);
      m_counts[5]->Fill(counts[9]);
      m_counts[6]->Fill(counts[10]);
      m_counts[7]->Fill(counts[11]);
      m_counts[8]->Fill(counts[12]);
      m_counts[9]->Fill(counts[13]);
      m_counts[10]->Fill(counts[1]+counts[4]+counts[6]+counts[7]+counts[8]+counts[9]+counts[10]+counts[11]+counts[12]+counts[13]);
      int total = 0;
      for(int i=0; i<64; i++)
         total = total + counts[i];
      m_counts[11]->Fill(total);
      
      // Fill anode relevant data
      // A1 - A2
      // [0] = A1B - A2B
      // [1] = A1B - A2M
      // [2] = A1B - A2T
      // [3] = A1M - A2B
      // [4] = A1M - A2M
      // [5] = A1M - A2T
      // [6] = A1T - A2B
      // [7] = A1T - A2M
      // [8] = A1T - A2T
      if (multiplicity > 0) {
         for (int i=0; i < 3; i++) {
            if(earliest_times[i] < 999990) {
               if(earliest_times[3] < 999990)
                  a1_a2_diff[0+3*i]->Fill(earliest_times[i]-earliest_times[3]);
               if(earliest_times[4] < 999990)
                  a1_a2_diff[1+3*i]->Fill(earliest_times[i]-earliest_times[4]);
               if(earliest_times[5] < 999990)
                  a1_a2_diff[2+3*i]->Fill(earliest_times[i]-earliest_times[5]);
            }
         }
      }
      
      // Fill cathode relevant data
      for(int i = 0; i < 4; i++) {
         int index = 6+i*2;
         // Only fill data for straight paths, multiplicity of 1
         //      if (apm == 1 || apm == 16 || apm == 256) {
         if (multiplicity > 0) {
            // First find a valid anode time for corresponding anode
            // This should be the smallest time between the 3 earliest_time,
            // since multiplicity = 1, the other 2 would be 999999
            double anode_time = 999999.0;
            for (int j=0; j<3; j++) {
               if (i == 0 || i == 1) {
                  if (earliest_times[j] < anode_time)
                     anode_time = earliest_times[j];
               }
               else {
                  if (earliest_times[3+j] < anode_time)
                     anode_time = earliest_times[3+j];
               }
            }
            // Check both cathode data are non-zero
            if(earliest_times[index] < 999999.0 && earliest_times[index+1] < 999999.0) {
               // Fill x/y_sum - 2*anode_time
               x_y_sum[i]->Fill(earliest_times[index+1]+earliest_times[index] - 2*anode_time);
               // R-L or B-T
               if(i==0 || i==2) {
                  x_y_diff[i]->Fill(earliest_times[index+1]-earliest_times[index]);
                  x_y_diff_vs_sum_silly[i]->Fill(earliest_times[index+1]+earliest_times[index] - 2*anode_time, earliest_times[index+1]-earliest_times[index]);
               }
               else {
                  x_y_diff[i]->Fill(earliest_times[index]-earliest_times[index+1]);
                  x_y_diff_vs_sum_silly[i]->Fill(earliest_times[index+1]+earliest_times[index] - 2*anode_time, earliest_times[index]-earliest_times[index+1]);
               }
            }
         }
      }
   } //end UpdateHistograms

   void BeginRun(TARunInfo* runinfo)
   {
      printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      fCounter = 0;
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      //fATX->BeginRun(runinfo->fRunNo);
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("EndRun, run %d, events %d\n", runinfo->fRunNo, fCounter);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
      //fATX->EndRun();
      //char fname[1024];
      //sprintf(fname, "output%05d.pdf", runinfo->fRunNo);
      //fATX->fH->fCanvas->SaveAs(fname);
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id != 1)
         return flow;

      TMBank* b = event->FindBank("EMMT");

      if (!b)
         return flow;

      int bklen = b->data_size;
      const char* bkptr = event->GetBankData(b);
         
      if (!bkptr)
         return flow;

      printf("EMMA TDC, pointer: %p, len %d\n", bkptr, bklen);

      while (bklen > 0) {
         v1190event *te = UnpackV1190(&bkptr, &bklen, fModule->fVerboseV1190);
         if (te == NULL)
            break;
         te->Print();
         UpdateHistograms(runinfo, te);
         delete te;
      }

      time_t now = time(NULL);

      static time_t t = 0;

      if (now - t > 15) {
         t = now;
         PlotHistograms(runinfo);
      }

      fCounter++;
      fModule->fTotalEventCounter++;

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

void EmmaModule::Init(const std::vector<std::string> &args)
{
   printf("Init!\n");

   fVerboseV1190 = false;

   for (unsigned i=0; i<args.size(); i++) {
      if (args[i] == "--verbose-v1190")
         fVerboseV1190 = true;
   }

   fTotalEventCounter = 0;
   TARootHelper::fgDir->cd(); // select correct ROOT directory
}
   
void EmmaModule::Finish()
{
   printf("Finish!\n");
   printf("Counted %d events grand total\n", fTotalEventCounter);
}
   
TARunInterface* EmmaModule::NewRun(TARunInfo* runinfo)
{
   printf("NewRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
   return new EmmaRun(runinfo, this);
}

TARegisterModule tarm(new EmmaModule);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
