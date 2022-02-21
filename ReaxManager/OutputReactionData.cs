using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
namespace ReaxManager
{
    public class OutputReactionData
    {
        private AtomInputData atomInputData;
        private List<int> totalMolNum;
        private List<int> totalSpecies;
        private List<List<int>> totalTargetMoleculeNum;
        private List<List<List<List<string>>>> totalReaction;
        List<Dictionary<string, int>> totalSpeciesDict;
        List<Dictionary<string, string>> totalSmilesPlusH;

        public OutputReactionData(AtomInputData atomInputData, List<int> totalMolNum, List<int> totalSpecies, List<List<int>> totalTargetMoleculeNum, 
                                  List<List<List<List<string>>>> totalReaction, List<Dictionary<string, int>> totalSpeciesDict, List<Dictionary<string, string>> totalSmilesPlusH)
        {
            this.atomInputData = atomInputData;
            this.totalMolNum = totalMolNum;
            this.totalSpecies = totalSpecies;
            this.totalTargetMoleculeNum = totalTargetMoleculeNum;
            this.totalReaction = totalReaction;
            this.totalSpeciesDict = totalSpeciesDict;
            this.totalSmilesPlusH = totalSmilesPlusH;
        }

        public void OutputToTextFile()
        {
            try
            {
                File.WriteAllText(@"ReactionData.txt", "Reaction Data" + Environment.NewLine);
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Molecules" + Environment.NewLine);
            for (int i = 0; i < atomInputData.TotalTimeStep - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalMolNum[i]}" + Environment.NewLine);
            }

            for(int i = 0; i < atomInputData.TargetMoleculeSmilesList.Count; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"Number of {atomInputData.TargetMoleculeSmilesList[0][i]}" + Environment.NewLine);
                for (int j = 0; j < atomInputData.TotalTimeStep - 1; j++)
                {
                    File.AppendAllText(@"ReactionData.txt", $"{totalTargetMoleculeNum[j][i]}" + Environment.NewLine);
                }
            }            

            File.AppendAllText(@"ReactionData.txt", "Number of Species" + Environment.NewLine);
            for (int i = 0; i < atomInputData.TotalTimeStep - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalSpecies[i]}" + Environment.NewLine);
            }
           
            File.AppendAllText(@"ReactionData.txt", "List of Species" + Environment.NewLine);
            for (int i = 0; i < atomInputData.TotalTimeStep - 1; i += atomInputData.TotalTimeStep-2)
            {                
                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                foreach (KeyValuePair<string, int> smiles in totalSpeciesDict[i])
                {
                    string molFormula = "";
                    int cCount, hCount, fCount, oCount;
                    cCount = 0; hCount = 0; fCount = 0; oCount = 0;
                    foreach (char atomName in totalSmilesPlusH[i][smiles.Key])
                    {
                        switch (atomName)
                        {
                            case 'C':
                                cCount++;
                                break;
                            case 'H':
                                hCount++;
                                break;
                            case 'F':
                                fCount++;
                                break;
                            case 'O':
                                oCount++;
                                break;                            
                        }
                    }                    
                    molFormula += (cCount != 0) ? ("C" + cCount.ToString()) : "";
                    molFormula += (hCount != 0) ? ("H" + hCount.ToString()) : "";
                    molFormula += (fCount != 0) ? ("F" + fCount.ToString()) : "";
                    molFormula += (oCount != 0) ? ("O" + oCount.ToString()) : "";
                    if (i == 0 || i == atomInputData.TotalTimeStep - 2)
                    {
                        File.AppendAllText(@"ReactionData.txt", $"●{smiles.Key} \"({molFormula})\": {smiles.Value}" + Environment.NewLine);
                    }
                }               
            }          

            File.AppendAllText(@"ReactionData.txt", "ReactionPerTime" + Environment.NewLine);
            for (int i = 0; i < atomInputData.TotalTimeStep - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                for (int i2 = 0; i2 < totalReaction[i].Count; i2++)
                {
                    File.AppendAllText(@"ReactionData.txt", $"------反応{i2 + 1}------" + Environment.NewLine);
                    for (int j = 0; j < totalReaction[i][i2].Count; j++)
                    {
                        for (int k = 0; k < totalReaction[i][i2][j].Count; k++)
                        {
                            File.AppendAllText(@"ReactionData.txt", $"{totalReaction[i][i2][j][k]}" + Environment.NewLine);
                            if (k == totalReaction[i][i2][j].Count - 1)
                            {
                                break;
                            }
                            File.AppendAllText(@"ReactionData.txt", "+" + Environment.NewLine);
                        }
                        if (j == 0)
                        {
                            File.AppendAllText(@"ReactionData.txt", "--↓↓--" + Environment.NewLine);
                        }
                    }
                    File.AppendAllText(@"ReactionData.txt", "------------" + Environment.NewLine);
                }
            }
        }
    }
}
