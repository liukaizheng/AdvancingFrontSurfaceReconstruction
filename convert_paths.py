import json
import os

with open('./compile_commands.json', 'r') as f:
    loadStr = json.load(f)

for i, data in enumerate(loadStr):
    command = data['command']
    cmdList = command.split(' ')
    change = False
    for j, cmdStr in enumerate(cmdList):
        if cmdStr:
            if cmdStr[0] == '@':
                change = True
                cmdStr = cmdStr[1:]
                includeFile = os.path.join(data['directory'], cmdStr)
                with open(includeFile, 'r') as incFile:
                    cmdStr = incFile.read()
                    cmdStr = cmdStr.replace('\n',  ' ')
                    cmdStr.rstrip()
                    cmdList[j] = cmdStr
    if change:
        command = ""
        for cmdStr in cmdList:
            if cmdStr:
                command += cmdStr + ' '
        command.rstrip()
        loadStr[i]['command'] = command

with open('./compile_commands.json', 'w') as f:
    json.dump(loadStr, f, indent = 4)





