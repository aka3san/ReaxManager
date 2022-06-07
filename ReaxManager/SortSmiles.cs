using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ReaxManager
{
    public class SortSmiles
    {
        private List<string> molNumTosmiles;
        private List<string> molNumToString;
        public SortSmiles(List<string> molNumToSmiles, List<string> molNumToString)
        {
            this.molNumTosmiles = new List<string>(molNumToSmiles);
            this.molNumToString = new List<string>(molNumToString);
        }

        public void SortSpecies(ref int totalSpecies, Dictionary<string, int> stringToSpeciesNum, Dictionary<string, string> smilesPlusH)
        {
            List<string> smiles = new List<string>(molNumTosmiles);
            List<int> smilesCount = new();
            for (int i = 0; i < smiles.Count; i++)
            {
                smiles[i] = smiles[i].Replace("[", "");
                smiles[i] = smiles[i].Replace("]", "");
                smiles[i] = smiles[i].Replace(":", "");
                smiles[i] = smiles[i].Replace("1", "");
                smiles[i] = smiles[i].Replace("H", "");
                if (smiles[i] == "=O")
                {
                    continue;
                }
                string smiles_temp = smiles[i];
                smiles[i] = smiles[i].Replace("=", "");
                smiles[i] = smiles[i].Replace("(", "");
                smiles[i] = smiles[i].Replace(")", "");
                if (!smilesCount.Contains(smiles[i].Count()))
                {
                    smilesCount.Add(smiles[i].Count());
                    stringToSpeciesNum.Add(smiles_temp, 1);
                    smilesPlusH.Add(smiles_temp, molNumToString[i]);
                }
                else
                {
                    foreach (string key in stringToSpeciesNum.Keys)
                    {
                        if (smiles[i].Length == key.Replace("=", "").Replace("(", "").Replace(")", "").Length)
                        {
                            stringToSpeciesNum[key] += 1;
                        }
                    }
                }
            }
            totalSpecies = stringToSpeciesNum.Count;
        }

        //それぞれの分子を特定して、数え上げていく
        public List<int> CountTargetMoleculeNumbers(List<string> targetMoleculeSmilesList)
        {
            List<string> smilesList = new List<string>(molNumTosmiles);
            List<int> molecularWeightList = GetMolecularWeightList(targetMoleculeSmilesList);
            List<int> moleculeNumbers = new List<int>();
            for (int i=0; i<molecularWeightList.Count; i++)
            {
                moleculeNumbers.Add(0);
            }
            
            foreach (string smiles in smilesList)
            {
                if (smiles.Length > 1000)
                {
                    continue;
                }

                string smiles2 = smiles.Replace("H", "");
                smiles2 = smiles2.Replace("[", "");
                smiles2 = smiles2.Replace("]", "");
                smiles2 = smiles2.Replace(":", "");
                smiles2 = smiles2.Replace("1", "");
                smiles2 = smiles2.Replace("(", "");
                smiles2 = smiles2.Replace(")", "");                
                smiles2 = smiles2.Replace("=", "");
                for(int i=0; i<molecularWeightList.Count; i++)
                {
                    //if(smiles2.Length % molecularWeightList[i] == 0 && smiles2.Length >= molecularWeightList[i])
                    if (smiles2.Length == molecularWeightList[i])
                    {
                        moleculeNumbers[i] = moleculeNumbers[i] + smiles2.Length/ molecularWeightList[i];
                    }
                }                
            }
            return moleculeNumbers;
        }

        private List<int> GetMolecularWeightList(List<string> smilesList)
        {
            List<int> molecularWeightList = new List<int>();
            List<string> smilesList_copy = new List<string>(smilesList);
            for(int i=0; i<smilesList_copy.Count; i++)
            {
                smilesList_copy[i] = smilesList_copy[i].Replace("H", "");
                smilesList_copy[i] = smilesList_copy[i].Replace("[", "");
                smilesList_copy[i] = smilesList_copy[i].Replace("]", "");
                smilesList_copy[i] = smilesList_copy[i].Replace(":", "");
                smilesList_copy[i] = smilesList_copy[i].Replace("1", "");
                smilesList_copy[i] = smilesList_copy[i].Replace("(", "");
                smilesList_copy[i] = smilesList_copy[i].Replace(")", "");
                smilesList_copy[i] = smilesList_copy[i].Replace("=", "");
                molecularWeightList.Add(smilesList_copy[i].Count());
            }
            return molecularWeightList;
        }
    }
}
