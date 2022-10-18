import sys
import os
import openpyxl
import xlwings as xw
import csv
import shutil

def pathStrip(input,levels=0): # cd .. from file path input for as many levels as specified
        output = input
        for i in range (levels + 1):
            output, tmp2 = os.path.split(output)
        return output


acceptanceCriteria = {1: {12: "Acceptance Criteria", 13: ""},
                      2: {12: "\"Constants\" column contains \"Constant 1\" for IV1_A samples only", 13: "=IF(COUNTIFS(E2:E151, \"Constant 1\", C2:C151, \"IV1_A\")=COUNTIFS(C2:C151, \"IV1_A\"), TRUE, FALSE)"},
                      3: {12: "\"Constants\" column contains \"Constant 3\" for IV1_B-IV2_A samples only", 13: "=IF(COUNTIFS(E2:E151, \"Constant 3\", C2:C151, \"IV1_B\", D2:D151, \"IV2_A\")=COUNTIFS(C2:C151, \"IV1_B\", D2:D151, \"IV2_A\"),TRUE, FALSE)"},
                      4: {12: "\"Constants\" column contain \"N/A\" for IV1_B-IV2_B samples", 13: "=IF(COUNTIFS(E2:E151, \"N/A\", C2:C151, \"IV1_B\", D2:D151, \"IV2_B\")=COUNTIFS(C2:C151, \"IV1_B\", D2:D151, \"IV2_B\"), TRUE, FALSE)"},
                      5: {12: "\"Constant 2\" column contains \"Constant Value\" for IV1_A and IV1_B-IV2_A samples", 13: "=IF(AND(COUNTIFS(F2:F151, \"Constant Value\", C2:C151, \"IV1_A\")=COUNTIFS(C2:C151, \"IV1_A\"), COUNTIFS(F2:F151, \"Constant Value\", C2:C151, \"IV1_B\", D2:D151, \"IV2_A\")=COUNTIFS(C2:C151, \"IV1_B\", D2:D151, \"IV2_A\")), TRUE, FALSE)"},
                      6: {12: "\"Constant 2\" column contains \"N/A\" for IV1_B-IV2_B samples", 13: "=IF(COUNTIFS(F2:F151,\"N/A\",C2:C151,\"IV1_B\",D2:D151,\"IV2_B\")=COUNTIFS(C2:C151,\"IV1_B\",D2:D151,\"IV2_B\"),TRUE,FALSE)"},
                      7: {12: "\"Rand Vars\" column contains an RV1 random variable for IV1_A-IV2_A samples only", 13: "=IF(COUNTIFS(G2:G151, \"*RV1*\", C2:C151, \"IV1_A\", D2:D151, \"IV2_A\")=COUNTIFS(C2:C151, \"IV1_A\", D2:D151, \"IV2_A\"), TRUE, FALSE)"},
                      8: {12: "\"Rand Vars\" column contains an RV3 random variable for IV1_B samples only", 13: "=IF(COUNTIFS(G2:G151, \"*RV3*\", C2:C151, \"IV1_B\")=COUNTIFS( C2:C151, \"IV1_B\"), TRUE, FALSE)"},
                      9: {12: "\"Rand Vars\" column contains \"N/A\" for IV1_A-IV2_B samples", 13: "=IF(COUNTIFS(G2:G151,\"N/A\",C2:C151,\"IV1_A\",D2:D151,\"IV2_B\")=COUNTIFS(C2:C151,\"IV1_A\",D2:D151,\"IV2_B\"),TRUE,FALSE)"},
                      10: {12: "\"Rand Var 2\" column contains a number between 1-5", 13: "=SUMPRODUCT(--(H2:H151={1,2,3,4,5}))=COUNT(H2:H151)"},
                      11: {12: "\"Only A-A\" column contains a blank for IV1_A-IV2_A samples only", 13: "=IF(COUNTIFS(I2:I151, \" \", C2:C151, \"IV1_A\", D2:D151, \"IV2_A\")=COUNTIFS(C2:C151, \"IV1_A\", D2:D151, \"IV2_A\"), TRUE, FALSE)"},
                      12: {12: "\"Only A-A\" column contains \"N/A\" for IV1_A-IV2_B and IV1_B samples", 13: "=IF(AND(COUNTIFS(I2:I151, \"N/A\", C2:C151, \"IV1_A\", D2:D151, \"IV2_B\")=COUNTIFS(C2:C151, \"IV1_A\", D2:D151, \"IV2_B\"),COUNTIFS(I2:I151, \"N/A\", C2:C151, \"IV1_B\")=COUNTIFS(C2:C151, \"IV1_B\")), TRUE, FALSE)"},
                      13: {12: "\"Only IV2_B\" column contains a blank for IV2_B samples only", 13: "=IF(COUNTIFS(J2:J151, \" \", D2:D151, \"IV2_B\")=COUNTIFS(D2:D151, \"IV2_B\"), TRUE, FALSE)"},
                      14: {12: "\"Only IV2_B\" column contains \"N/A\" for IV2_A samples only", 13: "=IF(COUNTIFS(J2:J151, \"N/A\", D2:D151, \"IV2_A\")=COUNTIFS( D2:D151, \"IV2_A\"), TRUE, FALSE)"},
                      15: {12: "", 13: "=IF(SUMIF(M2:M14, TRUE)=COUNT(M2:M14), \"Pass\", \"Fail\")"}}

cd = os.getcwd()
sampleSpaceFilePath = cd + "\\" + sys.argv[1]
variablesFilePath = cd + "\\" + sys.argv[2]
iterations = eval(sys.argv[3])

testResults = []
sampleTableFile = pathStrip(variablesFilePath)
if not os.path.exists(sampleTableFile + "\\ValidationCSVs"):
    os.mkdir(sampleTableFile + "\\ValidationCSVs")
if not os.path.exists(sampleTableFile + "\\ValidationWorkbooks"):
    os.mkdir(sampleTableFile + "\\ValidationWorkbooks")

for i in range(iterations):
    sampleTableName = "Val Sample Table " + str(i)
    os.system("python .\TestMatrixGenerator.py \"" + sampleSpaceFilePath + "\" \"" + variablesFilePath + "\" \"" + sampleTableName + "\"")
    sampleTableFilePath = sampleTableFile + "\\" + sampleTableName + ".csv"
    sampleTableCSVFilePath = sampleTableFile + "\\ValidationCSVs\\" + sampleTableName + ".csv"
    shutil.move(sampleTableFilePath, sampleTableCSVFilePath)
    #os.system("move \"" + sampleTableFilePath + "\" \"" + sampleTableCSVFilePath + "\"")
    sampleTableWorkbookFilePath = sampleTableFile + "\\ValidationWorkbooks\\Val Sample Table " + str(i) + ".xlsx"
    book = openpyxl.Workbook()
    sheet = book.active
    with open(sampleTableCSVFilePath) as f:
        reader = csv.reader(f)
        for row in reader:
            sheet.append(row)
    cells = sheet['L1':'M15']
    for row in range(1,16):
        for col in [12, 13]:
            sheet.cell(row, col).value = acceptanceCriteria[row][col]
    book.save(sampleTableWorkbookFilePath)
    book.close()
    app = xw.App(visible=False)
    book = app.books.open(sampleTableWorkbookFilePath)
    testResults.append(book.sheets['Sheet'].range('M15').value)
    book.save()
    book.close()
    app.quit()
    

print("Pass" if all(result == "Pass" for result in testResults) else "Fail")
