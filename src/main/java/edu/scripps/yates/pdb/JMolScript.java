package edu.scripps.yates.pdb;

import java.util.ArrayList;
import java.util.List;

public class JMolScript {
	private final List<String> commandList = new ArrayList<String>();

	public JMolScript() {

	}

	public JMolScript(String command) {
		addCommand(command);
	}

	public JMolScript addCommand(int command) {
		return addCommand(String.valueOf(command));
	}

	public JMolScript addCommand(String command) {
		if (command.endsWith(";")) {
			// remove the last ";"
			command = command.substring(0, command.length() - 2);
		}
		commandList.add(command);
		return this;
	}

	public String getCommandsToExecute() {
		StringBuilder sb = new StringBuilder();
		for (String command : commandList) {
			sb.append(command).append(";");
		}
		return sb.toString();
	}

	public String getCommandsToExecuteIndifferentLines() {
		StringBuilder sb = new StringBuilder();
		for (String command : commandList) {
			sb.append(command).append(";\n");
		}
		return sb.toString();
	}

	public void appendToLastCommand(String commandToAppend) {
		final int lastIndex = commandList.size() - 1;
		final String previousCommand = commandList.get(lastIndex);
		commandList.set(lastIndex, previousCommand + " " + commandToAppend);
	}
}
